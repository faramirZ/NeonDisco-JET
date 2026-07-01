#!/bin/bash
set -euo pipefail

#### JET Identification Pipeline — Master Orchestrator
#### Runs Step1 (STAR) -> Step2 (JET calling) -> Step3 (netMHCpan) -> Step4 (aggregation)
#### Supports multiple samples in a single run via a metadata file.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LIB_DIR="${SCRIPT_DIR}/lib"
source "${LIB_DIR}/common.sh"

print_banner() {
cat << "EOF"
       _ ______ ______   ____  _              ___
      | |  ____|__   __| |  _ \(_)            / _ \
      | | |__     | |    | |_) |_ _ __   ___ | | | |
  _   | |  __|    | |    |  _ <| | '_ \ / _ \| | | |
 | |__| | |____   | |    | |_) | | |_) |  __/| |_| |
  \____/|______|  |_|    |____/|_| .__/ \___| \___/
                                  | |
                                  |_|
EOF
echo -e "  JET Identification Pipeline — v1.0\n"
}

usage() {
cat << EOF
Usage: $(basename "$0") --samples FILE [OPTIONS]

Required:
  --samples FILE        Tab-separated metadata file. One row per sample:
                         name<TAB>R1.fastq.gz<TAB>R2.fastq.gz<TAB>HLA_alleles<TAB>repeatmasker_file
                         Example row:
                         379T   /data/379T_R1.fastq.gz   /data/379T_R2.fastq.gz   HLA-A11:01,HLA-A02:03   /meta/379T.repeatmasker_reformat.txt

Reference/Genome options:
  --genome GENOME        Genome build: hg19 | hg38 | mm9 | mm10 (default: hg38)
  --genome-fasta PATH    Genome FASTA (required unless --skip-step1)
  --gtf PATH             Gene GTF file (required unless --skip-step1 and --skip-step2)
  --star-index DIR       STAR genome index directory (required unless --skip-step1)

Tool paths:
  --star-bin DIR         Directory containing STAR binary
  --samtools-bin DIR     Directory containing samtools binary
  --rscript-bin PATH     Path to Rscript binary
  --r-script PATH        Path to JET_analysis_filtered.R
  --netmhcpan-bin DIR    Directory containing netMHCpan binary

Run control:
  --outputs-dir DIR      Output directory root (default: ./outputData)
  --tmp-dir DIR          Temp directory for STAR (default: ./tmpData)
  --threads N            Threads for STAR/samtools (default: 8)
  --read-length N         STAR sjdbOverhang (default: 100)
  --skip-step1           Skip STAR alignment (Step 1 output must already exist)
  --skip-step2           Skip JET identification (Step 2 output must already exist)
  --skip-step3           Skip netMHCpan (Step 3 output must already exist)
  --skip-step4           Skip binder aggregation
  -h, --help             Show this help

Example:
  $(basename "$0") --samples samples.tsv \\
      --genome hg38 \\
      --genome-fasta /data/ref/GRCh38.fa \\
      --gtf /data/ref/GRCh38.115.gtf \\
      --star-index /data/star_index \\
      --star-bin /opt/.pixi/envs/default/bin \\
      --samtools-bin /opt/.pixi/envs/default/bin \\
      --rscript-bin /opt/R/4.5.3/bin/Rscript \\
      --r-script ./JET_analysis_filtered.R \\
      --netmhcpan-bin /opt/netMHCpan-4.2/Linux_x86_64/bin
EOF
exit 1
}

#-------------------------------------------------------
##-- Default values
#-------------------------------------------------------

SAMPLES_FILE=""
GENOME="hg38"
GENOME_FASTA=""
GTF_FILE=""
STAR_INDEX_DIR=""
STAR_BIN=""
SAMTOOLS_BIN=""
RSCRIPT_BIN=""
R_SCRIPT=""
NETMHCPAN_BIN=""
OUTPUTS_DIR="$(pwd)/outputData"
TMP_DIR="$(pwd)/tmpData"
THREADS=8
READ_LENGTH=100
SKIP_STEP1=false
SKIP_STEP2=false
SKIP_STEP3=false
SKIP_STEP4=false

#-------------------------------------------------------
##-- Argument parsing
#-------------------------------------------------------

[ "$#" -eq 0 ] && usage

while [ "$#" -gt 0 ]; do
    case "$1" in
        --samples)        SAMPLES_FILE="$2"; shift 2 ;;
        --genome)          GENOME="$2"; shift 2 ;;
        --genome-fasta)    GENOME_FASTA="$2"; shift 2 ;;
        --gtf)             GTF_FILE="$2"; shift 2 ;;
        --star-index)      STAR_INDEX_DIR="$2"; shift 2 ;;
        --star-bin)        STAR_BIN="$2"; shift 2 ;;
        --samtools-bin)    SAMTOOLS_BIN="$2"; shift 2 ;;
        --rscript-bin)     RSCRIPT_BIN="$2"; shift 2 ;;
        --r-script)        R_SCRIPT="$2"; shift 2 ;;
        --netmhcpan-bin)   NETMHCPAN_BIN="$2"; shift 2 ;;
        --outputs-dir)     OUTPUTS_DIR="$2"; shift 2 ;;
        --tmp-dir)         TMP_DIR="$2"; shift 2 ;;
        --threads)         THREADS="$2"; shift 2 ;;
        --read-length)     READ_LENGTH="$2"; shift 2 ;;
        --skip-step1)      SKIP_STEP1=true; shift ;;
        --skip-step2)      SKIP_STEP2=true; shift ;;
        --skip-step3)      SKIP_STEP3=true; shift ;;
        --skip-step4)      SKIP_STEP4=true; shift ;;
        -h|--help)         usage ;;
        *) log_error "Unknown argument: $1"; usage ;;
    esac
done

print_banner

#-------------------------------------------------------
##-- Validate arguments
#-------------------------------------------------------

[ -n "${SAMPLES_FILE}" ] || die "Missing required argument: --samples"
require_file "${SAMPLES_FILE}" "Samples metadata file"

if [ "${SKIP_STEP1}" = false ]; then
    [ -n "${GENOME_FASTA}" ]  || die "Missing --genome-fasta (required unless --skip-step1)"
    [ -n "${GTF_FILE}" ]      || die "Missing --gtf (required unless --skip-step1)"
    [ -n "${STAR_INDEX_DIR}" ] || die "Missing --star-index (required unless --skip-step1)"
    [ -n "${STAR_BIN}" ]      || die "Missing --star-bin (required unless --skip-step1)"
    [ -n "${SAMTOOLS_BIN}" ]  || die "Missing --samtools-bin (required unless --skip-step1)"
    require_file "${GENOME_FASTA}" "Genome FASTA"
    require_file "${GTF_FILE}" "GTF file"
fi

if [ "${SKIP_STEP2}" = false ]; then
    [ -n "${GTF_FILE}" ]      || die "Missing --gtf (required unless --skip-step2)"
    [ -n "${RSCRIPT_BIN}" ]   || die "Missing --rscript-bin (required unless --skip-step2)"
    [ -n "${R_SCRIPT}" ]      || die "Missing --r-script (required unless --skip-step2)"
    require_file "${RSCRIPT_BIN}" "Rscript binary"
    require_file "${R_SCRIPT}" "R script"
fi

if [ "${SKIP_STEP3}" = false ]; then
    [ -n "${NETMHCPAN_BIN}" ] || die "Missing --netmhcpan-bin (required unless --skip-step3)"
fi

case "${GENOME}" in
    hg19|hg38|mm9|mm10) ;;
    *) die "Invalid --genome: ${GENOME} (must be hg19, hg38, mm9, or mm10)" ;;
esac

mkdir -p "${OUTPUTS_DIR}" "${TMP_DIR}"

export STAR_BIN SAMTOOLS_BIN GENOME_FASTA GTF_FILE STAR_INDEX_DIR OUTPUTS_DIR TMP_DIR
export READ_LENGTH THREADS RSCRIPT_BIN R_SCRIPT GENOME NETMHCPAN_BIN

LOG_DIR="${OUTPUTS_DIR}/logs"
mkdir -p "${LOG_DIR}"
MASTER_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

log_info "Logging full run to: ${MASTER_LOG}"
log_info "Samples file        : ${SAMPLES_FILE}"
log_info "Genome              : ${GENOME}"
log_info "Threads             : ${THREADS}"
log_info "Steps to skip       : step1=${SKIP_STEP1} step2=${SKIP_STEP2} step3=${SKIP_STEP3} step4=${SKIP_STEP4}"

#-------------------------------------------------------
##-- Validate samples file content up front
#-------------------------------------------------------

N_SAMPLES=0
while IFS=$'\t' read -r NAME R1 R2 HLA REPEATS; do
    [ -z "${NAME}" ] && continue
    [[ "${NAME}" =~ ^#.* ]] && continue
    N_SAMPLES=$((N_SAMPLES+1))
    if [ "${SKIP_STEP1}" = false ]; then
        require_file "${R1}" "R1 fastq for sample ${NAME}"
        require_file "${R2}" "R2 fastq for sample ${NAME}"
    fi
    if [ "${SKIP_STEP2}" = false ]; then
        require_file "${REPEATS}" "RepeatMasker file for sample ${NAME}"
    fi
    if [ "${SKIP_STEP3}" = false ]; then
        echo "${HLA}" | grep -qE "HLA-[ABC]" || die "Sample ${NAME}: malformed HLA_alleles field: '${HLA}'"
    fi
done < "${SAMPLES_FILE}"

[ "${N_SAMPLES}" -gt 0 ] || die "No valid sample rows found in ${SAMPLES_FILE}"
log_info "Found ${N_SAMPLES} sample(s) to process"

#-------------------------------------------------------
##-- Run pipeline per sample
#-------------------------------------------------------

FAILED_SAMPLES=()
SAMPLE_NUM=0

while IFS=$'\t' read -r NAME R1 R2 HLA REPEATS; do
    [ -z "${NAME}" ] && continue
    [[ "${NAME}" =~ ^#.* ]] && continue
    SAMPLE_NUM=$((SAMPLE_NUM+1))

    echo ""
    echo "================================================================"
    log_step "SAMPLE [${SAMPLE_NUM}/${N_SAMPLES}]: ${NAME}"
    echo "================================================================"

    SAMPLE_OK=true

    if [ "${SKIP_STEP1}" = false ]; then
        OUT_SAMPLE_DIR=$("${SCRIPT_DIR}/step1_alignment.sh" "${NAME}" "${R1}" "${R2}" | tail -1) \
            || { log_error "Step 1 failed for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step1"); SAMPLE_OK=false; }
    else
        OUT_SAMPLE_DIR=$(ls -td "${OUTPUTS_DIR}/${NAME}"_* 2>/dev/null | head -1)
        [ -n "${OUT_SAMPLE_DIR}" ] || { log_error "No existing output dir for ${NAME} (--skip-step1 set)"; FAILED_SAMPLES+=("${NAME}:step1-missing"); SAMPLE_OK=false; }
    fi

    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP2}" = false ]; then
        "${SCRIPT_DIR}/step2_jet_identification.sh" "${NAME}" "${OUT_SAMPLE_DIR}" "${REPEATS}" \
            || { log_error "Step 2 failed for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step2"); SAMPLE_OK=false; }
    fi

    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP3}" = false ]; then
        "${SCRIPT_DIR}/step3_netmhcpan.sh" "${NAME}" "${OUT_SAMPLE_DIR}" "${HLA}" \
            || { log_error "Step 3 failed for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step3"); SAMPLE_OK=false; }
    fi

    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP4}" = false ]; then
        "${SCRIPT_DIR}/aggregate_binders.sh" "${NAME}" "${OUT_SAMPLE_DIR}" \
            || { log_error "Step 4 failed for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step4"); SAMPLE_OK=false; }
    fi

    if [ "${SAMPLE_OK}" = true ]; then
        log_ok "Sample ${NAME} completed successfully. Output: ${OUT_SAMPLE_DIR}"
    fi

done < "${SAMPLES_FILE}" 2>&1 | tee -a "${MASTER_LOG}"

#-------------------------------------------------------
##-- Final summary
#-------------------------------------------------------

echo ""
echo "================================================================"
if [ "${#FAILED_SAMPLES[@]}" -eq 0 ]; then
    log_ok "Pipeline finished. All ${N_SAMPLES} sample(s) completed successfully."
else
    log_error "Pipeline finished with ${#FAILED_SAMPLES[@]} failure(s):"
    for f in "${FAILED_SAMPLES[@]}"; do
        log_error "  - ${f}"
    done
    log_error "See log for details: ${MASTER_LOG}"
    exit 1
fi
echo "================================================================"
