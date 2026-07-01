#!/bin/bash
set -euo pipefail

#### JET Pipeline — Step 3: netMHCpan binding prediction
#### Usage: step3_netmhcpan.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>
#### Requires env vars: NETMHCPAN_BIN, LIB_DIR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR> <HLA_ALLELES>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2
HLA_ALLELES=$3

: "${NETMHCPAN_BIN:?NETMHCPAN_BIN not set}"

require_dir  "${OUT_SAMPLE_DIR}" "Step 2 output directory"
require_file "${NETMHCPAN_BIN}/netMHCpan" "netMHCpan binary"
echo "${HLA_ALLELES}" | grep -qE "HLA-[ABC]" || die "HLA_ALLELES malformed: ${HLA_ALLELES}"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"

log_step "STEP 3 — netMHCpan: ${SAMPLE_NAME}"
log_info "HLA alleles: ${HLA_ALLELES}"

for SIZE in 8 9 10 11; do
    FASTA_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.fasta"
    if [ ! -f "${FASTA_FILE}" ]; then
        log_warn "FASTA not found for size ${SIZE}, skipping: ${FASTA_FILE}"
        continue
    fi

    OUT_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.netmhcpan4.txt"
    TMP_FILE="${OUT_FILE}.tmp"

    run_cmd "netMHCpan (${SAMPLE_NAME}, size=${SIZE})" \
        bash -c "'${NETMHCPAN_BIN}/netMHCpan' -a '${HLA_ALLELES}' -f '${FASTA_FILE}' -l '${SIZE}' -xls -xlsfile '${OUT_FILE}' -inptype 0 > '${TMP_FILE}'"
done

log_ok "Step 3 complete for ${SAMPLE_NAME} (sizes 8-11)."
