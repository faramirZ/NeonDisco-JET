#!/bin/bash
set -euo pipefail

#### NeonDisco-TEAL — Step 3: pVACbind MHC-I binding prediction
#### Usage: step3_pvacbind.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>
#### Requires env vars exported by master: PVACTOOLS_SIF, THREADS, LIB_DIR
####
#### HLA_ALLELES: comma-separated, accepts both formats:
####   pVACbind format  — HLA-A*11:01 (preferred)
####   netMHCpan format — HLA-A11:01  (converted internally)
####
#### pVACbind handles unknown/unsupported alleles ggez — it logs a
#### warning and skips them. If ALL alleles are invalid, this script fails.
####
#### Output per size: <OUT_SAMPLE_DIR>/pvacbind_size<N>/MHC_Class_I/
####   <SAMPLE>.MHC_I.all_epitopes.tsv          — all predictions
####   <SAMPLE>.MHC_I.filtered.tsv              — pVACbind default IC50 filter
####   <SAMPLE>.MHC_I.all_epitopes.aggregated.tsv
####
#### Optional env vars:
####   KEEP_TMP=true   — keep pVACbind tmp files after run (default: remove)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2
HLA_ALLELES=$3

: "${PVACTOOLS_SIF:?PVACTOOLS_SIF not set}"
: "${THREADS:?THREADS not set}"

KEEP_TMP="${KEEP_TMP:-false}"

require_dir  "${OUT_SAMPLE_DIR}" "Step 2 output directory"
require_file "${PVACTOOLS_SIF}"  "pVACtools SIF image"

echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]" \
    || die "HLA_ALLELES malformed — expected comma-separated HLA-A/B/C alleles, got: ${HLA_ALLELES}"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"

log_step "STEP 3 — pVACbind: ${SAMPLE_NAME}"
log_info "SIF image            : ${PVACTOOLS_SIF}"
log_info "HLA alleles (raw)    : ${HLA_ALLELES}"
log_info "Keep tmp files       : ${KEEP_TMP}"

# =============================================================================
#  CONVERT HLA ALLELE FORMAT IF NEEDED
#  samples.tsv may use netMHCpan format (HLA-A11:01)
#  pVACbind requires                    HLA-A*11:01 (asterisk after gene letter)
#  If asterisk already present, skip conversion.
# =============================================================================
if echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]\*"; then
    PVAC_HLA="${HLA_ALLELES}"
    log_info "HLA alleles (pVACbind): ${PVAC_HLA} (already correct, no conversion needed)"
else
    PVAC_HLA=$(echo "${HLA_ALLELES}" | sed 's/HLA-\([ABC]\)\([0-9]\)/HLA-\1*\2/g')
    log_info "HLA alleles (pVACbind): ${PVAC_HLA} (converted from netMHCpan format)"
fi

# =============================================================================
#  ALGORITHMS
#  All confirmed available in pVACtools.sif on crmy-penguin.
#  pVACbind handles unsupported alleles per-algorithm gracefully — it logs
#  "Allele X not valid. Skipping." and continues with remaining alleles.
#  No pre-validation needed.
# =============================================================================
ALGORITHMS="NetMHCpan NetMHCpanEL MHCflurryEL BigMHC_IM PRIME DeepImmuno"
log_info "Algorithms           : ${ALGORITHMS}"
echo ""

# =============================================================================
#  KEEP-TMP FLAG
#  -k is the correct short flag confirmed from pvacbind run --help
# =============================================================================
if [ "${KEEP_TMP}" = "true" ]; then
    TMP_FLAG="-k"
    log_info "Tmp files will be KEPT (KEEP_TMP=true)"
else
    TMP_FLAG=""
    log_info "Tmp files will be REMOVED after each size run (set KEEP_TMP=true to keep)"
fi
echo ""

# =============================================================================
#  RUN pVACbind PER PEPTIDE SIZE (8, 9, 10, 11)
#  One run per size so outputs map cleanly to per-size ids.txt for Step 4.
#
#  apptainer bind notes:
#    -B /home/faramir  — covers input FASTAs, output dirs, pipeline scripts
#    -B /GB            — covers data on /GB partition (symlinked paths)
#    --no-mount tmp    — container uses its own /tmp; pVACbind writes
#                        temp files to pvacbind_size${N}/MHC_Class_I/tmp/
#                        not to /tmp, so this is safe
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

    PVAC_OUT_DIR="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}"
    mkdir -p "${PVAC_OUT_DIR}"

    EXPECTED_OUT="${PVAC_OUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.all_epitopes.tsv"

    log_info "[size=${SIZE}] Input FASTA : ${FASTA_FILE}"
    log_info "[size=${SIZE}] Output dir  : ${PVAC_OUT_DIR}"
    log_info "[size=${SIZE}] Alleles     : ${PVAC_HLA}"

    run_cmd "pVACbind (${SAMPLE_NAME}, size=${SIZE})" \
        apptainer exec \
            --no-mount tmp \
            -B /GB:/GB \
            -B /home/faramir:/home/faramir \
            "${PVACTOOLS_SIF}" \
            pvacbind run \
                "${FASTA_FILE}" \
                "${SAMPLE_NAME}" \
                "${PVAC_HLA}" \
                ${ALGORITHMS} \
                "${PVAC_OUT_DIR}" \
                -e1 "${SIZE}" \
                -t "${THREADS}" \
                ${TMP_FLAG}

    # Confirm output landed and has content
    if [ -f "${EXPECTED_OUT}" ]; then
        NLINES=$(tail -n +2 "${EXPECTED_OUT}" | wc -l)
        if [ "${NLINES}" -eq 0 ]; then
            log_warn "[size=${SIZE}] pVACbind ran but produced 0 predictions — all alleles may have been invalid or no peptides passed"
            SIZES_SKIPPED=$(( SIZES_SKIPPED + 1 ))
        else
            log_ok "[size=${SIZE}] pVACbind complete — ${NLINES} epitope predictions written"
            SIZES_RUN=$(( SIZES_RUN + 1 ))
        fi
    else
        log_warn "[size=${SIZE}] Expected output not found: ${EXPECTED_OUT}"
        log_warn "[size=${SIZE}] All HLA alleles may have been invalid for this run"
        SIZES_SKIPPED=$(( SIZES_SKIPPED + 1 ))
    fi
    echo ""
done

# =============================================================================
#  FINAL CHECK AND SUMMARY
# =============================================================================
if [ "${SIZES_RUN}" -eq 0 ]; then
    die "pVACbind produced no output for any size for ${SAMPLE_NAME} — check HLA alleles and FASTA files"
fi

log_ok "Step 3 complete for ${SAMPLE_NAME}."
log_info "  Sizes run     : ${SIZES_RUN}"
log_info "  Sizes skipped : ${SIZES_SKIPPED}"
log_info "  Alleles used  : ${PVAC_HLA}"
log_info "  Algorithms    : ${ALGORITHMS}"
