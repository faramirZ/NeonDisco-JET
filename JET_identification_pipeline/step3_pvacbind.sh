#!/bin/bash
set -euo pipefail

#### NeonDisco-TEAL — Step 3: pVACbind MHC-I binding prediction
#### Usage: step3_pvacbind.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>
#### Requires env vars exported by master: PVACTOOLS_SIF, THREADS, LIB_DIR
####
#### HLA_ALLELES: comma-separated in samples.tsv format (HLA-A11:01)
####              converted internally to pVACbind format  (HLA-A*11:01)
####
#### Output per size: <OUT_SAMPLE_DIR>/pvacbind_size<N>/MHC_Class_I/
####   <SAMPLE>.MHC_I.all_epitopes.tsv   — all predictions
####   <SAMPLE>.MHC_I.filtered.tsv       — pVACbind default IC50 filter
####   <SAMPLE>.MHC_I.all_epitopes.aggregated.tsv

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2
HLA_ALLELES=$3

: "${PVACTOOLS_SIF:?PVACTOOLS_SIF not set}"
: "${THREADS:?THREADS not set}"

require_dir  "${OUT_SAMPLE_DIR}" "Step 2 output directory"
require_file "${PVACTOOLS_SIF}"  "pVACtools SIF image"

echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]" \
    || die "HLA_ALLELES malformed — expected comma-separated HLA-A/B/C alleles, got: ${HLA_ALLELES}"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"

log_step "STEP 3 — pVACbind: ${SAMPLE_NAME}"
log_info "SIF image            : ${PVACTOOLS_SIF}"
log_info "HLA alleles (raw)    : ${HLA_ALLELES}"

# =============================================================================
#  CONVERT HLA ALLELE FORMAT
#  samples.tsv uses netMHCpan format:  HLA-A11:01
#  pVACbind requires:                  HLA-A*11:01 (asterisk after gene letter)
# =============================================================================
# Convert HLA format if needed: HLA-A11:01 → HLA-A*11:01
# Checks first — if asterisk already present, skips conversion entirely
if echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]\*"; then
    PVAC_HLA="${HLA_ALLELES}"
    log_info "HLA alleles (pVACbind format): ${PVAC_HLA} (already correct, no conversion needed)"
else
    PVAC_HLA=$(echo "${HLA_ALLELES}" | sed 's/HLA-\([ABC]\)\([0-9]\)/HLA-\1*\2/g')
    log_info "HLA alleles (pVACbind format): ${PVAC_HLA} (converted from netMHCpan format)"
fi

# =============================================================================
#  ALGORITHMS
#  All confirmed available in pVACtools.sif image v4.0.10 (2023-12-20) on crmy-penguin
# =============================================================================
ALGORITHMS="NetMHCpan NetMHCpanEL MHCflurryEL BigMHC_IM PRIME DeepImmuno"
log_info "Algorithms           : ${ALGORITHMS}"
echo ""

# =============================================================================
#  RUN pVACbind PER PEPTIDE SIZE (8, 9, 10, 11)
#  One run per size so outputs map cleanly to per-size ids.txt for Step 4
# =============================================================================
for SIZE in 8 9 10 11; do
    FASTA_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.fasta"
    if [ ! -f "${FASTA_FILE}" ]; then
        log_warn "FASTA not found for size ${SIZE}, skipping: ${FASTA_FILE}"
        continue
    fi

    PVAC_OUT_DIR="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}"
    mkdir -p "${PVAC_OUT_DIR}"

    EXPECTED_OUT="${PVAC_OUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.all_epitopes.tsv"

    log_info "[size=${SIZE}] Input FASTA : ${FASTA_FILE}"
    log_info "[size=${SIZE}] Output dir  : ${PVAC_OUT_DIR}"

    # apptainer exec notes:
    # --no-mount tmp  — container gets its own /tmp; pVACbind writes
    #                   temp files under ${PVAC_OUT_DIR}/MHC_Class_I/tmp/
    #                   (not /tmp) so this is safe
    # -B /home/faramir — binds all input FASTAs, output dirs, and the
    #                   pipeline scripts in one mount since everything
    #                   lives under /home/faramir on crmy-penguin.
    #                   If OUTPUTS_DIR or inputData ever moves to /GB,
    #                   add: -B /GB:/GB
    run_cmd "pVACbind (${SAMPLE_NAME}, size=${SIZE})" \
        apptainer exec \
            --no-mount tmp \
            -B /home/faramir:/home/faramir \
            "${PVACTOOLS_SIF}" \
            pvacbind run \
                "${FASTA_FILE}" \
                "${SAMPLE_NAME}" \
                "${PVAC_HLA}" \
                ${ALGORITHMS} \
                "${PVAC_OUT_DIR}" \
                -e1 "${SIZE}" \
                -t "${THREADS}"

    # Confirm expected output landed
    if [ -f "${EXPECTED_OUT}" ]; then
        NLINES=$(tail -n +2 "${EXPECTED_OUT}" | wc -l)
        log_ok "[size=${SIZE}] pVACbind complete — ${NLINES} epitope predictions written"
    else
        log_warn "[size=${SIZE}] Expected output not found: ${EXPECTED_OUT}"
    fi
    echo ""
done

log_ok "Step 3 complete for ${SAMPLE_NAME} (sizes 8-11)."
log_ok "Alleles used: ${PVAC_HLA}"
