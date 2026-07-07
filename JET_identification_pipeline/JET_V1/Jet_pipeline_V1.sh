#!/bin/bash
# =============================================================================
#  neondisco-jet — JET Identification Pipeline
#  Master Orchestrator
#
#  Created by Alexandre Houy and Christel Goudot
#  Modified and adapted by Muhammad Faramir
#
#  Usage:
#    bash run_neondisco_jet.sh --samples <samples.txt> --hla <alleles> [OPTIONS]
#
#  Reference:
#    Burbage & Rocañín-Arjó et al., Science Immunology 2023
# =============================================================================

set -euo pipefail

# =============================================================================
#  COLOUR HELPERS
# =============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

log()     { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${BOLD}$*${RESET}"             | tee -a "${MASTER_LOG}"; }
log_ok()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${GREEN}✔  $*${RESET}"         | tee -a "${MASTER_LOG}"; }
log_err() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${RED}✘  ERROR: $*${RESET}"    | tee -a "${MASTER_LOG}"; }
log_warn(){ echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] ${YELLOW}⚠  WARNING: $*${RESET}" | tee -a "${MASTER_LOG}"; }
log_step(){ echo -e "\n${CYAN}${BOLD}════════════════════════════════════════${RESET}"  | tee -a "${MASTER_LOG}"
            echo -e "${CYAN}${BOLD}  $*${RESET}"                                        | tee -a "${MASTER_LOG}"
            echo -e "${CYAN}${BOLD}════════════════════════════════════════${RESET}\n"  | tee -a "${MASTER_LOG}"; }

# =============================================================================
#  BANNER
# =============================================================================
print_banner() {
    echo -e "${BOLD}"
    echo "  ███╗   ██╗███████╗ ██████╗ ███╗   ██╗██████╗ ██╗███████╗ ██████╗  ██╗     ██╗███████╗████████╗"
    echo "  ████╗  ██║██╔════╝██╔═══██╗████╗  ██║██╔══██╗██║██╔════╝██╔════╝ ██  ██╗ ██║██╔════╝╚══██╔══╝"
    echo "  ██╔██╗ ██║█████╗  ██║   ██║██╔██╗ ██║██║  ██║██║███████╗██║      ██  ██║ ██║█████╗     ██║   "
    echo "  ██║╚██╗██║██╔══╝  ██║   ██║██║╚██╗██║██║  ██║██║╚════██║██║      ██  ██║ ██║██╔══╝     ██║   "
    echo "  ██║ ╚████║███████╗╚██████╔╝██║ ╚████║██████╔╝██║███████║╚██████╗ ╚████╔╝ ██║███████╗   ██║   "
    echo "  ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚═════╝ ╚═╝╚══════╝ ╚═════╝  ╚═══╝  ╚═╝╚══════╝   ╚═╝   "
    echo -e "${RESET}"
    echo -e "  ${CYAN}JET Identification Pipeline — neondisco-jet${RESET}"
    echo -e "  ${CYAN}Burbage & Rocañín-Arjó et al., Science Immunology 2023${RESET}"
    echo -e "  ${CYAN}Adapted for human GRCh38 by Muhammad Faramir${RESET}\n"
}

# =============================================================================
#  USAGE
# =============================================================================
usage() {
    echo -e "${BOLD}Usage:${RESET}"
    echo "  bash $(basename "$0") --samples <samples.txt> --hla <alleles> [OPTIONS]"
    echo ""
    echo -e "${BOLD}Required arguments:${RESET}"
    echo "  --samples   Path to samples.txt (3 columns: R1_path R2_path sample_name)"
    echo "  --hla       Comma-separated HLA alleles in netMHCpan format"
    echo "              e.g. HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04"
    echo ""
    echo -e "${BOLD}Optional arguments:${RESET}"
    echo "  --steps     Which steps to run: 1, 2, 3, or any combination e.g. '1,2' or '2,3'"
    echo "              Default: 1,2,3 (run all steps)"
    echo "  --skip-step1  Skip Step 1 (use if STAR alignment already done)"
    echo "  --help      Show this message"
    echo ""
    echo -e "${BOLD}Example (full pipeline):${RESET}"
    echo "  bash $(basename "$0") \\"
    echo "    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.txt \\"
    echo "    --hla \"HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04\""
    echo ""
    echo -e "${BOLD}Example (skip Step 1, run Steps 2 and 3 only):${RESET}"
    echo "  bash $(basename "$0") \\"
    echo "    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.txt \\"
    echo "    --hla \"HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04\" \\"
    echo "    --steps 2,3"
    echo ""
    exit 1
}

# =============================================================================
#  ARGUMENT PARSING
# =============================================================================
SAMPLES_FILE=""
HLA_ALLELES=""
RUN_STEPS="1,2,3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)   SAMPLES_FILE="$2";  shift 2 ;;
        --hla)       HLA_ALLELES="$2";   shift 2 ;;
        --steps)     RUN_STEPS="$2";     shift 2 ;;
        --help|-h)   print_banner; usage ;;
        *)           echo "Unknown argument: $1"; usage ;;
    esac
done

# =============================================================================
#  VALIDATE REQUIRED ARGUMENTS
# =============================================================================
print_banner

if [ -z "${SAMPLES_FILE}" ]; then
    echo -e "${RED}Error: --samples is required.${RESET}"; usage
fi
if [ -z "${HLA_ALLELES}" ]; then
    echo -e "${RED}Error: --hla is required.${RESET}"; usage
fi
if [ ! -f "${SAMPLES_FILE}" ]; then
    echo -e "${RED}Error: samples file not found: ${SAMPLES_FILE}${RESET}"; exit 1
fi
if ! echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]"; then
    echo -e "${RED}Error: HLA alleles do not appear valid. Expected format: HLA-A11:01,HLA-B13:01${RESET}"; exit 1
fi

# =============================================================================
#  CONFIGURATION — EDIT THESE PATHS FOR YOUR ENVIRONMENT
# =============================================================================

# Software paths
STAR_BIN=/home/faramir/repos/neondisco-jet/.pixi/envs/default/bin
SAMTOOLS_BIN=/home/faramir/repos/neondisco-jet/.pixi/envs/default/bin
R_BIN=/opt/R/4.5.3/bin
NETMHCPAN_DIR=/home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2

# Pipeline script directory
PIPELINE_DIR=/home/faramir/repos/neondisco-jet/JET_identification_pipeline
R_SCRIPT=${PIPELINE_DIR}/Copy_JET_analysis_filtered.R

# Reference files
FASTA_FILE=/home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF_FILE=/home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf

# Output and log directories
OUTPUT_DIR=/home/faramir/data/neondisco-jet/outputData
LOG_DIR=/home/faramir/data/neondisco-jet/logDir
TMP_DIR=/home/faramir/data/neondisco-jet/ErrorDir
STAR_INDEX_DIR=/home/faramir/data/neondisco-jet/starIndexesDir

# Pipeline settings
GENOME="GRCh38"
READ_LENGTH=100
THREADS=16
PEPTIDE_SIZES="8 9 10 11"

# =============================================================================
#  SETUP LOGGING
# =============================================================================
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}" "${TMP_DIR}"
RUN_TIMESTAMP=$(date +"%Y%m%d_%Hh%Mm%Ss")
MASTER_LOG="${LOG_DIR}/pipeline_run_${RUN_TIMESTAMP}.log"
touch "${MASTER_LOG}"

# =============================================================================
#  VALIDATE TOOLS AND REFERENCE FILES
# =============================================================================
log_step "PRE-FLIGHT CHECKS"

check_tool() {
    local tool_path="$1"
    local tool_name="$2"
    if [ ! -f "${tool_path}" ] && ! command -v "${tool_path}" &>/dev/null; then
        log_err "${tool_name} not found at: ${tool_path}"
        exit 1
    fi
    log_ok "${tool_name} found"
}

check_file() {
    local file_path="$1"
    local file_desc="$2"
    if [ ! -f "${file_path}" ]; then
        log_err "${file_desc} not found: ${file_path}"
        exit 1
    fi
    log_ok "${file_desc} found"
}

check_tool "${STAR_BIN}/STAR"           "STAR aligner"
check_tool "${SAMTOOLS_BIN}/samtools"   "samtools"
check_tool "${R_BIN}/Rscript"           "Rscript"
check_tool "${NETMHCPAN_DIR}/netMHCpan" "netMHCpan"
check_file "${FASTA_FILE}"              "Reference FASTA"
check_file "${GTF_FILE}"                "Reference GTF"
check_file "${R_SCRIPT}"                "JET R script"
check_file "${SAMPLES_FILE}"            "Samples file"

log ""
log "Pipeline log       : ${MASTER_LOG}"
log "Samples file       : ${SAMPLES_FILE}"
log "HLA alleles        : ${HLA_ALLELES}"
log "Steps to run       : ${RUN_STEPS}"
log "Peptide sizes      : ${PEPTIDE_SIZES}"
log "Output directory   : ${OUTPUT_DIR}"

# =============================================================================
#  HELPER FUNCTION: ELAPSED TIME
# =============================================================================
elapsed() {
    local start=$1
    local end=$(date +%s)
    local secs=$(( end - start ))
    printf "%02d:%02d:%02d" $(( secs/3600 )) $(( (secs%3600)/60 )) $(( secs%60 ))
}

# =============================================================================
#  STEP 1 — STAR ALIGNMENT
# =============================================================================
run_step1() {
    local name=$1
    local fastqR1=$2
    local fastqR2=$3

    log_step "STEP 1 — STAR ALIGNMENT | Sample: ${name}"
    local t_start=$(date +%s)

    # Derive uncompressed paths
    local fastqR1File="${fastqR1%.gz}"
    local fastqR2File="${fastqR2%.gz}"

    # Create output directory with today's date
    local run_day=$(date +"%Y%m%d")
    local sample_dir="${OUTPUT_DIR}/${name}_${run_day}"
    mkdir -p "${sample_dir}"

    local prefix="${sample_dir}/${name}"

    # Decompress if needed
    if [ ! -f "${fastqR1File}" ] && [ -f "${fastqR1}" ]; then
        log "Decompressing ${fastqR1}..."
        gunzip -c "${fastqR1}" > "${fastqR1File}"
    fi
    if [ ! -f "${fastqR2File}" ] && [ -f "${fastqR2}" ]; then
        log "Decompressing ${fastqR2}..."
        gunzip -c "${fastqR2}" > "${fastqR2File}"
    fi

    # Validate fastq files exist
    if [ ! -f "${fastqR1File}" ] && [ ! -f "${fastqR1}" ]; then
        log_err "R1 fastq not found: ${fastqR1}"
        return 1
    fi
    if [ ! -f "${fastqR2File}" ] && [ ! -f "${fastqR2}" ]; then
        log_err "R2 fastq not found: ${fastqR2}"
        return 1
    fi

    log "Input R1        : ${fastqR1File}"
    log "Input R2        : ${fastqR2File}"
    log "Output prefix   : ${prefix}"
    log "STAR index      : ${STAR_INDEX_DIR}"

    # Run STAR alignment
    log "Running STAR alignment..."
    ${STAR_BIN}/STAR \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --runThreadN ${THREADS} \
        --genomeDir ${STAR_INDEX_DIR} \
        --sjdbGTFfile ${GTF_FILE} \
        --sjdbOverhang ${READ_LENGTH} \
        --readFilesIn ${fastqR1File} ${fastqR2File} \
        --outFileNamePrefix ${prefix}_ \
        --outTmpDir ${TMP_DIR}/STAR_${name}_${RUN_TIMESTAMP} \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --bamRemoveDuplicatesType UniqueIdentical \
        --outFilterMismatchNoverLmax 0.04 \
        --outMultimapperOrder Random \
        --outFilterMultimapNmax 1000 \
        --winAnchorMultimapNmax 1000 \
        --chimOutType WithinBAM \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        2>&1 | tee -a "${MASTER_LOG}"

    if [ $? -ne 0 ]; then
        log_err "STAR alignment failed for sample ${name}"
        return 1
    fi
    log_ok "STAR alignment complete"

    # Sort and index chimeric BAM
    log "Sorting and indexing chimeric BAM..."
    local bamChimeric="${prefix}_Chimeric.out.bam"
    local bamChimericSort="${prefix}_Chimeric.out.sort.bam"
    local bamMain="${prefix}_Aligned.sortedByCoord.out.bam"

    ${SAMTOOLS_BIN}/samtools sort -@ ${THREADS} -o "${bamChimericSort}" -O bam "${bamChimeric}" 2>&1 | tee -a "${MASTER_LOG}"
    ${SAMTOOLS_BIN}/samtools index "${bamChimericSort}" 2>&1 | tee -a "${MASTER_LOG}"
    ${SAMTOOLS_BIN}/samtools index "${bamMain}" 2>&1 | tee -a "${MASTER_LOG}"

    log_ok "BAM sorting and indexing complete"
    log_ok "Step 1 finished for sample ${name} in $(elapsed ${t_start})"

    # Return the output directory so Step 2 can find it
    echo "${sample_dir}"
}

# =============================================================================
#  STEP 2 — JET IDENTIFICATION (R SCRIPT)
# =============================================================================
run_step2() {
    local name=$1
    local sample_dir=$2

    log_step "STEP 2 — JET IDENTIFICATION | Sample: ${name}"
    local t_start=$(date +%s)

    local prefix="${sample_dir}/${name}"
    local chimericFile="${prefix}_Chimeric.out.junction"
    local junctionFile="${prefix}_SJ.out.tab"
    local logFinalOut="${prefix}_Log.final.out"

    # Validate STAR output files exist
    for f in "${chimericFile}" "${junctionFile}" "${logFinalOut}"; do
        if [ ! -f "${f}" ]; then
            log_err "Required file not found: ${f}"
            log_err "Make sure Step 1 ran successfully before Step 2."
            return 1
        fi
    done

    # Extract library size from STAR log
    local libsize
    libsize="$(grep "Uniquely mapped reads number" "${logFinalOut}" | sed -r 's/[[:space:]]+//g' | cut -d '|' -f2)"

    if [ -z "${libsize}" ]; then
        log_err "Could not extract library size from ${logFinalOut}"
        return 1
    fi

    log "Chimeric file   : ${chimericFile}"
    log "Junction file   : ${junctionFile}"
    log "Library size    : ${libsize}"
    log "Peptide sizes   : ${PEPTIDE_SIZES}"

    # Run R script for each peptide size
    local size_failed=0
    for size in ${PEPTIDE_SIZES}; do
        log ""
        log "[$(date '+%Y-%m-%d %H:%M:%S')] Running JET analysis — peptide size ${size}..."

        ${R_BIN}/Rscript "${R_SCRIPT}" \
            --chimeric "${chimericFile}" \
            --junction "${junctionFile}" \
            --genome "${GENOME}" \
            --size "${size}" \
            --libsize "${libsize}" \
            --prefix "${prefix}" \
            --verbose \
            2>&1 | tee -a "${MASTER_LOG}"

        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log_ok "JET analysis size ${size}: SUCCESS"
        else
            log_err "JET analysis size ${size}: FAILED"
            size_failed=1
        fi
    done

    if [ ${size_failed} -eq 1 ]; then
        log_warn "One or more sizes failed in Step 2 for sample ${name}"
        return 1
    fi

    log_ok "Step 2 finished for sample ${name} in $(elapsed ${t_start})"
}

# =============================================================================
#  STEP 3 — NETMHCPAN BINDING PREDICTION
# =============================================================================
run_step3() {
    local name=$1
    local sample_dir=$2

    log_step "STEP 3 — MHC-I BINDING PREDICTION | Sample: ${name}"
    local t_start=$(date +%s)

    local prefix="${sample_dir}/${name}"

    log "HLA alleles     : ${HLA_ALLELES}"
    log "Peptide sizes   : ${PEPTIDE_SIZES}"
    log "Output prefix   : ${prefix}"

    local size_failed=0
    for size in ${PEPTIDE_SIZES}; do
        log ""
        log "[$(date '+%Y-%m-%d %H:%M:%S')] Running netMHCpan — peptide size ${size}..."

        local fastaFile="${prefix}_Fusions_chim2.junc2.size${size}.fasta"
        local netmhcpanFile="${prefix}_Fusions_chim2.junc2.size${size}.netmhcpan4.txt"
        local netmhcpanTmp="${netmhcpanFile}.tmp"

        if [ ! -f "${fastaFile}" ]; then
            log_warn "FASTA not found for size ${size}: ${fastaFile} — skipping"
            continue
        fi

        ${NETMHCPAN_DIR}/netMHCpan \
            -a "${HLA_ALLELES}" \
            -f "${fastaFile}" \
            -l "${size}" \
            -xls \
            -xlsfile "${netmhcpanFile}" \
            -inptype 0 \
            > "${netmhcpanTmp}" 2>&1

        if [ $? -eq 0 ]; then
            log_ok "netMHCpan size ${size}: SUCCESS → ${netmhcpanFile}"
        else
            log_err "netMHCpan size ${size}: FAILED"
            size_failed=1
        fi
    done

    if [ ${size_failed} -eq 1 ]; then
        log_warn "One or more sizes failed in Step 3 for sample ${name}"
        return 1
    fi

    log_ok "Step 3 finished for sample ${name} in $(elapsed ${t_start})"
}

# =============================================================================
#  MAIN LOOP — PROCESS ALL SAMPLES
# =============================================================================
log_step "STARTING PIPELINE RUN"
log "Timestamp     : ${RUN_TIMESTAMP}"
log "Samples file  : ${SAMPLES_FILE}"
log "HLA alleles   : ${HLA_ALLELES}"
log "Steps to run  : ${RUN_STEPS}"

PIPELINE_START=$(date +%s)
SAMPLE_COUNT=0
SAMPLE_FAILED=0

while IFS=$'\t' read -r fastqR1 fastqR2 name || [[ -n "${name}" ]]; do
    # Skip empty lines and comment lines
    [[ -z "${name}" || "${name}" == \#* ]] && continue

    SAMPLE_COUNT=$(( SAMPLE_COUNT + 1 ))
    log_step "PROCESSING SAMPLE ${SAMPLE_COUNT}: ${name}"

    sample_dir=""
    sample_ok=true

    # ── STEP 1 ──────────────────────────────────────────────────
    if echo "${RUN_STEPS}" | grep -q "1"; then
        sample_dir=$(run_step1 "${name}" "${fastqR1}" "${fastqR2}") || {
            log_err "Step 1 failed for sample ${name} — skipping Steps 2 and 3"
            SAMPLE_FAILED=$(( SAMPLE_FAILED + 1 ))
            sample_ok=false
        }
    else
        # Find existing output directory if Step 1 is skipped
        sample_dir=$(ls -td "${OUTPUT_DIR}/${name}_"* 2>/dev/null | head -1)
        if [ -z "${sample_dir}" ]; then
            log_err "No existing output directory found for ${name} in ${OUTPUT_DIR}"
            log_err "Cannot skip Step 1 when no prior run exists."
            SAMPLE_FAILED=$(( SAMPLE_FAILED + 1 ))
            sample_ok=false
        else
            log_ok "Using existing output directory: ${sample_dir}"
        fi
    fi

    # ── STEP 2 ──────────────────────────────────────────────────
    if ${sample_ok} && echo "${RUN_STEPS}" | grep -q "2"; then
        run_step2 "${name}" "${sample_dir}" || {
            log_err "Step 2 failed for sample ${name} — skipping Step 3"
            SAMPLE_FAILED=$(( SAMPLE_FAILED + 1 ))
            sample_ok=false
        }
    fi

    # ── STEP 3 ──────────────────────────────────────────────────
    if ${sample_ok} && echo "${RUN_STEPS}" | grep -q "3"; then
        run_step3 "${name}" "${sample_dir}" || {
            log_err "Step 3 failed for sample ${name}"
            SAMPLE_FAILED=$(( SAMPLE_FAILED + 1 ))
            sample_ok=false
        }
    fi

    if ${sample_ok}; then
        log_ok "All requested steps completed for sample ${name}"
        log "Output files in: ${sample_dir}"
    fi

done < "${SAMPLES_FILE}"

# =============================================================================
#  FINAL SUMMARY
# =============================================================================
log_step "PIPELINE COMPLETE"
log "Total samples processed : ${SAMPLE_COUNT}"
log "Samples succeeded       : $(( SAMPLE_COUNT - SAMPLE_FAILED ))"
log "Samples failed          : ${SAMPLE_FAILED}"
log "Total runtime           : $(elapsed ${PIPELINE_START})"
log "Full log saved to       : ${MASTER_LOG}"

if [ ${SAMPLE_FAILED} -gt 0 ]; then
    log_err "${SAMPLE_FAILED} sample(s) had errors. Review log: ${MASTER_LOG}"
    exit 1
else
    log_ok "Pipeline finished successfully for all ${SAMPLE_COUNT} sample(s)."
fi
