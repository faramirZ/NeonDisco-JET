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
log_info "HLA alleles provided : ${HLA_ALLELES}"

# =============================================================================
#  VALIDATE HLA ALLELES AGAINST netMHCpan LIBRARY
#  Runs -listMHC once, checks each allele, drops any not found
# =============================================================================
log_info "Querying netMHCpan library for supported alleles..."

# Declare with safe default first to avoid set -u unbound errors.
# Capture both stdout and stderr since netMHCpan -listMHC may exit non-zero
# even on success — || true prevents set -e from killing the script.
MHC_LIST=""
MHC_LIST=$(cd "${NETMHCPAN_BIN}" && ./netMHCpan -listMHC 2>&1 | grep -E "^HLA-|^H-2") || true

if [ -z "${MHC_LIST}" ]; then
    die "netMHCpan -listMHC returned no alleles — check that ${NETMHCPAN_BIN}/netMHCpan can find its data files"
fi

MHC_LIST_COUNT=$(echo "${MHC_LIST}" | wc -l)
log_info "netMHCpan library loaded: ${MHC_LIST_COUNT} alleles available"

# Declare arrays with safe defaults before the loop
VALID_ALLELES=()
DROPPED_ALLELES=()

IFS=',' read -ra ALLELE_ARRAY <<< "${HLA_ALLELES}"
for ALLELE in "${ALLELE_ARRAY[@]}"; do
    # Strip all whitespace variants including carriage returns
    ALLELE=$(echo "${ALLELE}" | tr -d ' \t\r\n')
    [ -z "${ALLELE}" ] && continue

    # Use awk exact line matching — more reliable than grep -xF
    if echo "${MHC_LIST}" | awk -v a="${ALLELE}" '$0==a{found=1} END{exit (found!=1)}' 2>/dev/null; then
        VALID_ALLELES+=("${ALLELE}")
        log_ok  "  ✔  Found in library   : ${ALLELE}"
    else
        DROPPED_ALLELES+=("${ALLELE}")
        log_warn "  ✘  Not in library, dropping: ${ALLELE}"
    fi
done

# -- Summary of validation
echo ""
log_info "── HLA Allele Validation Summary for ${SAMPLE_NAME} ──"
log_info "  Provided  : ${#ALLELE_ARRAY[@]} allele(s)"
log_info "  Valid     : ${#VALID_ALLELES[@]} allele(s)"
log_info "  Dropped   : ${#DROPPED_ALLELES[@]} allele(s)"

if [ "${#DROPPED_ALLELES[@]}" -gt 0 ]; then
    log_warn "Dropped alleles (not in netMHCpan library):"
    for A in "${DROPPED_ALLELES[@]}"; do
        log_warn "    ✘  ${A}"
    done
fi

if [ "${#VALID_ALLELES[@]}" -eq 0 ]; then
    die "No valid HLA alleles remaining for ${SAMPLE_NAME} after filtering. Cannot run netMHCpan."
fi

# Build cleaned comma-separated allele string for netMHCpan
CLEAN_HLA=$(IFS=','; echo "${VALID_ALLELES[*]}")
log_ok "Proceeding with ${#VALID_ALLELES[@]} allele(s): ${CLEAN_HLA}"
echo ""

# =============================================================================
#  RUN netMHCpan PER PEPTIDE SIZE (8, 9, 10, 11)
# =============================================================================
for SIZE in 8 9 10 11; do
    FASTA_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.fasta"
    if [ ! -f "${FASTA_FILE}" ]; then
        log_warn "FASTA not found for size ${SIZE}, skipping: ${FASTA_FILE}"
        continue
    fi

    OUT_FILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.netmhcpan4.txt"
    TMP_FILE="${OUT_FILE}.tmp"

    run_cmd "netMHCpan (${SAMPLE_NAME}, size=${SIZE})" \
        bash -c "cd '${NETMHCPAN_BIN}' && ./netMHCpan -a '${CLEAN_HLA}' -f '${FASTA_FILE}' -l '${SIZE}' -xls -xlsfile '${OUT_FILE}' -inptype 0 > '${TMP_FILE}'"
done

log_ok "Step 3 complete for ${SAMPLE_NAME}."
log_ok "Alleles used: ${CLEAN_HLA}"
