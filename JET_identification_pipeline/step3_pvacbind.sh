#!/bin/bash
set -euo pipefail

#### JET Pipeline — Step 3: pVACbind MHC-I binding prediction
#### Replaces: step3_netmhcpan.sh
#### Usage: step3_pvacbind.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>
#### Requires env vars: PVACTOOLS_SIF, LIB_DIR
####
#### HLA allele format for pVACbind: HLA-A*11:01 (asterisk required)
#### samples.tsv must use this format — different from netMHCpan (HLA-A11:01)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2
HLA_ALLELES=$3

: "${PVACTOOLS_SIF:?PVACTOOLS_SIF not set — export path to pVACtools.sif}"

require_dir  "${OUT_SAMPLE_DIR}" "Step 2 output directory"
require_file "${PVACTOOLS_SIF}"  "pVACtools Apptainer image"

# Validate HLA format — pVACbind requires asterisk (HLA-A*11:01 not HLA-A11:01)
if ! echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]\*"; then
    die "HLA alleles must use pVACbind format with asterisk: HLA-A*11:01,HLA-B*13:01 — got: ${HLA_ALLELES}"
fi

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"
PVAC_OUT_DIR="${OUT_SAMPLE_DIR}/pvacbind_${SAMPLE_NAME}"

# Algorithms to run — two independent families, affinity + eluted ligand each
ALGORITHMS="NetMHCpan NetMHCpanEL MHCflurry MHCflurryEL"

log_step "STEP 3 — pVACbind: ${SAMPLE_NAME}"
log_info "HLA alleles  : ${HLA_ALLELES}"
log_info "Algorithms   : ${ALGORITHMS}"
log_info "Output dir   : ${PVAC_OUT_DIR}"
log_info "SIF image    : ${PVACTOOLS_SIF}"

# =============================================================================
#  VALIDATE HLA ALLELES AGAINST pVACbind
#  Check each allele is recognised before running the full prediction
# =============================================================================
log_info "Validating HLA alleles against pVACbind library..."

VALID_ALLELES=()
DROPPED_ALLELES=()

# Get supported alleles — strip \r from output (Windows line endings)
MHC_LIST_FILE=$(mktemp)
trap 'rm -f "${MHC_LIST_FILE}"' EXIT

apptainer exec \
    --no-mount tmp \
    -B /GB:/GB \
    -B /home/faramir:/home/faramir \
    "${PVACTOOLS_SIF}" \
    pvacbind valid_alleles 2>/dev/null \
    | tr -d '\r' > "${MHC_LIST_FILE}" \
    || die "Failed to query pVACbind valid alleles from ${PVACTOOLS_SIF}"

IFS=',' read -ra ALLELE_ARRAY <<< "${HLA_ALLELES}"
for ALLELE in "${ALLELE_ARRAY[@]}"; do
    ALLELE=$(echo "${ALLELE}" | tr -d ' \r')
    if grep -qxF "${ALLELE}" "${MHC_LIST_FILE}"; then
        VALID_ALLELES+=("${ALLELE}")
        log_ok  "  ✔  Found in library   : ${ALLELE}"
    else
        DROPPED_ALLELES+=("${ALLELE}")
        log_warn "  ✘  Not in library, dropping: ${ALLELE}"
    fi
done

echo ""
log_info "── HLA Allele Validation Summary for ${SAMPLE_NAME} ──"
log_info "  Provided  : ${#ALLELE_ARRAY[@]} allele(s)"
log_info "  Valid     : ${#VALID_ALLELES[@]} allele(s)"
log_info "  Dropped   : ${#DROPPED_ALLELES[@]} allele(s)"

if [ "${#DROPPED_ALLELES[@]}" -gt 0 ]; then
    log_warn "Dropped alleles (not in pVACbind library):"
    for A in "${DROPPED_ALLELES[@]}"; do
        log_warn "    ✘  ${A}"
    done
fi

if [ "${#VALID_ALLELES[@]}" -eq 0 ]; then
    die "No valid HLA alleles remaining for ${SAMPLE_NAME} after filtering. Cannot run pVACbind."
fi

# Build comma-separated allele string for pVACbind
CLEAN_HLA=$(IFS=','; echo "${VALID_ALLELES[*]}")
log_ok "Proceeding with ${#VALID_ALLELES[@]} allele(s): ${CLEAN_HLA}"
echo ""

# =============================================================================
#  RUN pVACbind PER PEPTIDE SIZE (8, 9, 10, 11)
#  pVACbind accepts all sizes in one call via --epitope-length
#  but we run per-size to match the Step 2 FASTA file structure
#  and keep output naming consistent with the old netMHCpan approach
# =============================================================================
SIZES_RUN=0
SIZES_SKIPPED=0

for SIZE in 8 9 10 11; do
    FASTA_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.fasta"

    if [ ! -f "${FASTA_FILE}" ]; then
        log_warn "FASTA not found for size ${SIZE}, skipping: ${FASTA_FILE}"
        SIZES_SKIPPED=$(( SIZES_SKIPPED + 1 ))
        continue
    fi

    SIZE_OUT_DIR="${PVAC_OUT_DIR}/size${SIZE}"
    mkdir -p "${SIZE_OUT_DIR}"

    log_info "Running pVACbind (${SAMPLE_NAME}, size=${SIZE})..."
    log_info "  Input FASTA : ${FASTA_FILE}"
    log_info "  Output dir  : ${SIZE_OUT_DIR}"

    run_cmd "pVACbind (${SAMPLE_NAME}, size=${SIZE})" \
        apptainer exec \
            --no-mount tmp \
            -B /GB:/GB \
            -B /home/faramir:/home/faramir \
            "${PVACTOOLS_SIF}" \
            pvacbind run \
                "${FASTA_FILE}" \
                "${SAMPLE_NAME}_size${SIZE}" \
                "${CLEAN_HLA}" \
                ${ALGORITHMS} \
                "${SIZE_OUT_DIR}" \
                --epitope-length "${SIZE}" \
                --binding-threshold 500 \
                --percentile-threshold 2 \
                --top-score-metric lowest \
                --keep-tmp-files

    # Verify output was produced
    ALL_EPITOPES="${SIZE_OUT_DIR}/MHC_Class_I/${SAMPLE_NAME}_size${SIZE}.all_epitopes.tsv"
    FILTERED="${SIZE_OUT_DIR}/MHC_Class_I/${SAMPLE_NAME}_size${SIZE}.filtered.tsv"

    if [ -f "${ALL_EPITOPES}" ]; then
        N_TOTAL=$(  tail -n +2 "${ALL_EPITOPES}" | wc -l)
        N_FILTERED=$(tail -n +2 "${FILTERED}"    | wc -l)
        log_ok "pVACbind (${SAMPLE_NAME}, size=${SIZE}) — ${N_TOTAL} epitopes, ${N_FILTERED} passed filter"
    else
        log_warn "pVACbind output not found for size ${SIZE}: ${ALL_EPITOPES}"
    fi

    SIZES_RUN=$(( SIZES_RUN + 1 ))
    echo ""
done

# =============================================================================
#  SUMMARY
# =============================================================================
log_ok "Step 3 complete for ${SAMPLE_NAME}."
log_info "  Sizes run     : ${SIZES_RUN}"
log_info "  Sizes skipped : ${SIZES_SKIPPED}"
log_info "  Alleles used  : ${CLEAN_HLA}"
log_info "  Algorithms    : ${ALGORITHMS}"
log_info "  Output root   : ${PVAC_OUT_DIR}"
