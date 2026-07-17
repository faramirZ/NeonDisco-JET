#!/bin/bash
set -euo pipefail

#### JET Pipeline — Step 2: JET identification (R-based)
#### Usage: step2_jet_identification.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <REPEATS_FILE>
#### Requires env vars: RSCRIPT_BIN, R_SCRIPT, GENOME, LIB_DIR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <REPEATS_FILE>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2
REPEATS_FILE=$3

: "${RSCRIPT_BIN:?RSCRIPT_BIN not set}"
: "${R_SCRIPT:?R_SCRIPT not set}"
: "${GENOME:?GENOME not set}"
: "${GTF_FILE:?GTF_FILE not set}"

require_dir  "${OUT_SAMPLE_DIR}" "Step 1 output directory"
require_file "${R_SCRIPT}" "R script"
require_file "${REPEATS_FILE}" "RepeatMasker file"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"
CHIMERIC_FILE="${PREFIX}_Chimeric.out.junction"
JUNCTION_FILE="${PREFIX}_SJ.out.tab"
LOG_FINAL_OUT="${PREFIX}_Log.final.out"

require_file "${CHIMERIC_FILE}" "Chimeric junction file"
require_file "${JUNCTION_FILE}" "SJ junction file"
require_file "${LOG_FINAL_OUT}" "STAR Log.final.out"

LIBSIZE=$(grep "Uniquely mapped reads number" "${LOG_FINAL_OUT}" | sed -r 's/[[:space:]]+//g' | cut -d '|' -f2)
[ -n "${LIBSIZE}" ] || die "Could not extract library size from ${LOG_FINAL_OUT}"

log_step "STEP 2 — JET identification: ${SAMPLE_NAME}"
log_info "Library size : ${LIBSIZE}"
log_info "Repeats file : ${REPEATS_FILE}"

for SIZE in 8 9 10 11; do
    run_cmd "R analysis (${SAMPLE_NAME}, size=${SIZE})" \
        "${RSCRIPT_BIN}" "${R_SCRIPT}" \
            --chimeric "${CHIMERIC_FILE}" \
            --junction "${JUNCTION_FILE}" \
            --genome "${GENOME}" \
            --size "${SIZE}" \
            --libsize "${LIBSIZE}" \
            --repeats "${REPEATS_FILE}" \
            --gtf "${GTF_FILE}" \
            --prefix "${PREFIX}" \
            --verbose
done

log_ok "Step 2 complete for ${SAMPLE_NAME} (sizes 8-11)."
