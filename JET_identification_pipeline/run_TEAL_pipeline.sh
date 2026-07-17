#!/bin/bash
set -euo pipefail

# =============================================================================
#  neondisco-teal — JET Identification Pipeline
#  Master Orchestrator
#
#  Created by Alexandre Houy and Christel Goudot
#  Modified and adapted by Muhammad Faramir
#
#  Usage:
#    bash run_neondisco_teal.sh --samples FILE [OPTIONS]
#
#  Reference:
#    Burbage & Rocañín-Arjó et al., Science Immunology 2023
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LIB_DIR="${SCRIPT_DIR}/lib"
source "${LIB_DIR}/common.sh"

# =============================================================================
#  COLOUR HELPERS
# =============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

log_info()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${BOLD}$*${RESET}"               | tee -a "${MASTER_LOG}"; }
log_ok()    { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${GREEN}✔  $*${RESET}"           | tee -a "${MASTER_LOG}"; }
log_error() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${RED}✘  ERROR: $*${RESET}"      | tee -a "${MASTER_LOG}"; }
log_warn()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${YELLOW}⚠  WARNING: $*${RESET}" | tee -a "${MASTER_LOG}"; }
log_step()  { echo -e "\n${CYAN}${BOLD}════════════════════════════════════════${RESET}" | tee -a "${MASTER_LOG}"
              echo -e "${CYAN}${BOLD}  $*${RESET}"                                       | tee -a "${MASTER_LOG}"
              echo -e "${CYAN}${BOLD}════════════════════════════════════════${RESET}\n" | tee -a "${MASTER_LOG}"; }

# =============================================================================
#  BANNER
# =============================================================================
print_banner() {
cat << "EOF"

  ███╗   ██╗███████╗ ██████╗ ███╗   ██╗██████╗ ██╗███████╗ ██████╗ ██████╗
  ████╗  ██║██╔════╝██╔═══██╗████╗  ██║██╔══██╗██║██╔════╝██╔════╝██╔═══██╗
  ██╔██╗ ██║█████╗  ██║   ██║██╔██╗ ██║██║  ██║██║███████╗██║     ██║   ██║
  ██║╚██╗██║██╔══╝  ██║   ██║██║╚██╗██║██║  ██║██║╚════██║██║     ██║   ██║
  ██║ ╚████║███████╗╚██████╔╝██║ ╚████║██████╔╝██║███████║╚██████╗╚██████╔╝
  ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚═════╝ ╚═╝╚══════╝ ╚═════╝ ╚═════╝

  ████████╗███████╗ █████╗ ██╗
  ╚══██╔══╝██╔════╝██╔══██╗██║
     ██║   █████╗  ███████║██║
     ██║   ██╔══╝  ██╔══██║██║
     ██║   ███████╗██║  ██║███████╗
     ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

EOF
    echo -e "  \033[0;36m\033[1mJET Identification Pipeline — neondisco-teal\033[0m"
    echo -e "  \033[0;36mBurbage & Rocañín-Arjó et al., Science Immunology 2023\033[0m"
    echo -e "  \033[0;36mAdapted for human GRCh38 by Muhammad Faramir\033[0m\n"
}

# =============================================================================
#  USAGE
# =============================================================================
usage() {
cat << EOF
Usage: $(basename "$0") --samples FILE [OPTIONS]

Required:
  --samples FILE        Tab-separated metadata file. One row per sample:
                         name<TAB>R1.fastq.gz<TAB>R2.fastq.gz<TAB>HLA_alleles<TAB>repeatmasker_file
                         Example row:
                         379T   /data/379T_R1.fastq   /data/379T_R2.fastq   HLA-A11:01,HLA-A02:03   /meta/379T.repeatmasker_reformat.txt

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
  --read-length N        STAR sjdbOverhang (default: 100)
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

# =============================================================================
#  DEFAULT VALUES
# =============================================================================
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

# =============================================================================
#  ARGUMENT PARSING
# =============================================================================
[ "$#" -eq 0 ] && usage

while [ "$#" -gt 0 ]; do
    case "$1" in
        --samples)        SAMPLES_FILE="$2"; shift 2 ;;
        --genome)         GENOME="$2"; shift 2 ;;
        --genome-fasta)   GENOME_FASTA="$2"; shift 2 ;;
        --gtf)            GTF_FILE="$2"; shift 2 ;;
        --star-index)     STAR_INDEX_DIR="$2"; shift 2 ;;
        --star-bin)       STAR_BIN="$2"; shift 2 ;;
        --samtools-bin)   SAMTOOLS_BIN="$2"; shift 2 ;;
        --rscript-bin)    RSCRIPT_BIN="$2"; shift 2 ;;
        --r-script)       R_SCRIPT="$2"; shift 2 ;;
        --netmhcpan-bin)  NETMHCPAN_BIN="$2"; shift 2 ;;
        --outputs-dir)    OUTPUTS_DIR="$2"; shift 2 ;;
        --tmp-dir)        TMP_DIR="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --read-length)    READ_LENGTH="$2"; shift 2 ;;
        --skip-step1)     SKIP_STEP1=true; shift ;;
        --skip-step2)     SKIP_STEP2=true; shift ;;
        --skip-step3)     SKIP_STEP3=true; shift ;;
        --skip-step4)     SKIP_STEP4=true; shift ;;
        -h|--help)        usage ;;
        *) echo "ERROR: Unknown argument: $1"; usage ;;
    esac
done

print_banner

# =============================================================================
#  SETUP LOG
#  (must come before any log_* calls that reference MASTER_LOG)
# =============================================================================
LOG_DIR="${OUTPUTS_DIR}/logs"
mkdir -p "${LOG_DIR}" "${OUTPUTS_DIR}" "${TMP_DIR}"
MASTER_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# =============================================================================
#  VALIDATE ARGUMENTS
# =============================================================================
[ -n "${SAMPLES_FILE}" ] || die "Missing required argument: --samples"
require_file "${SAMPLES_FILE}" "Samples metadata file"

if [ "${SKIP_STEP1}" = false ]; then
    [ -n "${GENOME_FASTA}" ]   || die "Missing --genome-fasta (required unless --skip-step1)"
    [ -n "${GTF_FILE}" ]       || die "Missing --gtf (required unless --skip-step1)"
    [ -n "${STAR_INDEX_DIR}" ] || die "Missing --star-index (required unless --skip-step1)"
    [ -n "${STAR_BIN}" ]       || die "Missing --star-bin (required unless --skip-step1)"
    [ -n "${SAMTOOLS_BIN}" ]   || die "Missing --samtools-bin (required unless --skip-step1)"
    require_file "${GENOME_FASTA}" "Genome FASTA"
    require_file "${GTF_FILE}" "GTF file"
fi

if [ "${SKIP_STEP2}" = false ]; then
    [ -n "${GTF_FILE}" ]    || die "Missing --gtf (required unless --skip-step2)"
    [ -n "${RSCRIPT_BIN}" ] || die "Missing --rscript-bin (required unless --skip-step2)"
    [ -n "${R_SCRIPT}" ]    || die "Missing --r-script (required unless --skip-step2)"
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

export STAR_BIN SAMTOOLS_BIN GENOME_FASTA GTF_FILE STAR_INDEX_DIR OUTPUTS_DIR TMP_DIR
export READ_LENGTH THREADS RSCRIPT_BIN R_SCRIPT GENOME NETMHCPAN_BIN

# =============================================================================
#  PIPELINE CONFIGURATION SUMMARY
# =============================================================================
log_info "════════════════════════════════════════"
log_info "  neondisco-teal — Pipeline Configuration"
log_info "════════════════════════════════════════"
log_info "Log file            : ${MASTER_LOG}"
log_info "Samples file        : ${SAMPLES_FILE}"
log_info "Genome              : ${GENOME}"
log_info "Threads             : ${THREADS}"
log_info "Read length         : ${READ_LENGTH}"
log_info "Outputs dir         : ${OUTPUTS_DIR}"
log_info "Tmp dir             : ${TMP_DIR}"
log_info "Skip Step 1 (STAR)  : ${SKIP_STEP1}"
log_info "Skip Step 2 (JET)   : ${SKIP_STEP2}"
log_info "Skip Step 3 (MHC)   : ${SKIP_STEP3}"
log_info "Skip Step 4 (Agg)   : ${SKIP_STEP4}"
log_info "════════════════════════════════════════"

# =============================================================================
#  PRE-FLIGHT: VALIDATE SAMPLES FILE CONTENT
# =============================================================================
log_step "PRE-FLIGHT: VALIDATING SAMPLES FILE"

N_SAMPLES=0
while IFS=$'\t' read -r NAME R1 R2 HLA REPEATS || [ -n "${NAME}" ]; do
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
    log_ok "Sample validated: ${NAME}"
done < "${SAMPLES_FILE}"

[ "${N_SAMPLES}" -gt 0 ] || die "No valid sample rows found in ${SAMPLES_FILE}"
log_info "Found ${N_SAMPLES} sample(s) to process"

# =============================================================================
#  RUN PIPELINE PER SAMPLE
# =============================================================================
FAILED_SAMPLES=()
SAMPLE_NUM=0
PIPELINE_START=$(date +%s)

while IFS=$'\t' read -r NAME R1 R2 HLA REPEATS; do
    [ -z "${NAME}" ] && continue
    [[ "${NAME}" =~ ^#.* ]] && continue
    SAMPLE_NUM=$((SAMPLE_NUM+1))

    log_step "SAMPLE [${SAMPLE_NUM}/${N_SAMPLES}]: ${NAME}"
    SAMPLE_START=$(date +%s)
    SAMPLE_OK=true

    # ── STEP 1: STAR ALIGNMENT ─────────────────────────────────────────────
    if [ "${SKIP_STEP1}" = false ]; then
        log_info "[Step 1/4] ── Starting STAR alignment for: ${NAME}"
        log_info "[Step 1/4]    R1 : ${R1}"
        log_info "[Step 1/4]    R2 : ${R2}"
        log_info "[Step 1/4]    Running STAR 2-pass alignment (this may take a while)..."

        STEP1_START=$(date +%s)
        OUT_SAMPLE_DIR=$("${SCRIPT_DIR}/step1_alignment.sh" "${NAME}" "${R1}" "${R2}" | tail -1) \
            || { log_error "Step 1 — STAR alignment FAILED for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step1"); SAMPLE_OK=false; }

        if [ "${SAMPLE_OK}" = true ]; then
            STEP1_END=$(date +%s)
            log_ok "[Step 1/4]    STAR alignment complete"
            log_ok "[Step 1/4]    Chimeric BAM sorting complete"
            log_ok "[Step 1/4]    BAM index files created"
            log_ok "[Step 1/4] ── Step 1 finished for ${NAME} ($(( STEP1_END - STEP1_START ))s) → ${OUT_SAMPLE_DIR}"
        fi
    else
        OUT_SAMPLE_DIR=$(ls -td "${OUTPUTS_DIR}/${NAME}"_* 2>/dev/null | head -1)
        if [ -n "${OUT_SAMPLE_DIR}" ]; then
            log_ok "[Step 1/4] Skipped — using existing output: ${OUT_SAMPLE_DIR}"
        else
            log_error "[Step 1/4] Skipped but no existing output directory found for ${NAME}"
            FAILED_SAMPLES+=("${NAME}:step1-missing")
            SAMPLE_OK=false
        fi
    fi

    # ── STEP 2: JET IDENTIFICATION ─────────────────────────────────────────
    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP2}" = false ]; then
        log_info "[Step 2/4] ── Starting JET identification for: ${NAME}"
        log_info "[Step 2/4]    RepeatMasker : ${REPEATS}"
        log_info "[Step 2/4]    Running R analysis for peptide sizes 8, 9, 10, 11..."

        STEP2_START=$(date +%s)
        "${SCRIPT_DIR}/step2_jet_identification.sh" "${NAME}" "${OUT_SAMPLE_DIR}" "${REPEATS}" \
            || { log_error "Step 2 — JET identification FAILED for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step2"); SAMPLE_OK=false; }

        if [ "${SAMPLE_OK}" = true ]; then
            STEP2_END=$(date +%s)
            log_ok "[Step 2/4]    Junction annotation complete"
            log_ok "[Step 2/4]    Fusion peptide FASTA files generated (sizes 8-11)"
            log_ok "[Step 2/4] ── Step 2 finished for ${NAME} ($(( STEP2_END - STEP2_START ))s)"
        fi
    elif [ "${SKIP_STEP2}" = true ]; then
        log_warn "[Step 2/4] Skipped for ${NAME}"
    fi

    # ── STEP 3: NETMHCPAN ──────────────────────────────────────────────────
    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP3}" = false ]; then
        log_info "[Step 3/4] ── Starting netMHCpan binding prediction for: ${NAME}"
        log_info "[Step 3/4]    HLA alleles : ${HLA}"
        log_info "[Step 3/4]    Running binding prediction for peptide sizes 8, 9, 10, 11..."

        STEP3_START=$(date +%s)
        "${SCRIPT_DIR}/step3_netmhcpan.sh" "${NAME}" "${OUT_SAMPLE_DIR}" "${HLA}" \
            || { log_error "Step 3 — netMHCpan FAILED for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step3"); SAMPLE_OK=false; }

        if [ "${SAMPLE_OK}" = true ]; then
            STEP3_END=$(date +%s)
            log_ok "[Step 3/4]    MHC-I binding prediction complete (sizes 8-11)"
            log_ok "[Step 3/4] ── Step 3 finished for ${NAME} ($(( STEP3_END - STEP3_START ))s)"
        fi
    elif [ "${SKIP_STEP3}" = true ]; then
        log_warn "[Step 3/4] Skipped for ${NAME}"
    fi

    # ── STEP 4: BINDER AGGREGATION ─────────────────────────────────────────
    if [ "${SAMPLE_OK}" = true ] && [ "${SKIP_STEP4}" = false ]; then
        log_info "[Step 4/4] ── Aggregating binders across all sizes for: ${NAME}"
        log_info "[Step 4/4]    Applying rank<2 (binder) and rank<0.5 (strong binder) filters..."

        STEP4_START=$(date +%s)
        "${SCRIPT_DIR}/aggregate_binders.sh" "${NAME}" "${OUT_SAMPLE_DIR}" \
            || { log_error "Step 4 — Aggregation FAILED for ${NAME}"; FAILED_SAMPLES+=("${NAME}:step4"); SAMPLE_OK=false; }

        if [ "${SAMPLE_OK}" = true ]; then
            STEP4_END=$(date +%s)
            log_ok "[Step 4/4]    Binder table written: ${OUT_SAMPLE_DIR}/all_sizes_binders.tsv"
            log_ok "[Step 4/4] ── Step 4 finished for ${NAME} ($(( STEP4_END - STEP4_START ))s)"
        fi
    elif [ "${SKIP_STEP4}" = true ]; then
        log_warn "[Step 4/4] Skipped for ${NAME}"
    fi

    # ── SAMPLE SUMMARY ─────────────────────────────────────────────────────
    SAMPLE_END=$(date +%s)
    SAMPLE_ELAPSED=$(printf "%02d:%02d:%02d" \
        $(( (SAMPLE_END-SAMPLE_START)/3600 )) \
        $(( ((SAMPLE_END-SAMPLE_START)%3600)/60 )) \
        $(( (SAMPLE_END-SAMPLE_START)%60 )))

    if [ "${SAMPLE_OK}" = true ]; then
        log_ok "══ Sample ${NAME} — all steps completed in ${SAMPLE_ELAPSED}"
        log_info "   Output directory: ${OUT_SAMPLE_DIR}"
    else
        log_error "══ Sample ${NAME} — finished with errors after ${SAMPLE_ELAPSED}"
        log_error "   Check log for details: ${MASTER_LOG}"
    fi

done < "${SAMPLES_FILE}" 2>&1 | tee -a "${MASTER_LOG}"

# =============================================================================
#  FINAL PIPELINE SUMMARY
# =============================================================================
PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
ELAPSED_FMT=$(printf "%02d:%02d:%02d" \
    $(( ELAPSED/3600 )) \
    $(( (ELAPSED%3600)/60 )) \
    $(( ELAPSED%60 )))

log_step "PIPELINE SUMMARY"
log_info "Total samples  : ${N_SAMPLES}"
log_info "Succeeded      : $(( N_SAMPLES - ${#FAILED_SAMPLES[@]} ))"
log_info "Failed         : ${#FAILED_SAMPLES[@]}"
log_info "Total runtime  : ${ELAPSED_FMT}"
log_info "Full log       : ${MASTER_LOG}"

if [ "${#FAILED_SAMPLES[@]}" -eq 0 ]; then
    log_ok "Pipeline finished — all ${N_SAMPLES} sample(s) completed successfully."
else
    log_error "Pipeline finished with ${#FAILED_SAMPLES[@]} failure(s):"
    for f in "${FAILED_SAMPLES[@]}"; do
        log_error "  ✘  ${f}"
    done
    log_error "Review the log for details: ${MASTER_LOG}"
    exit 1
fi
