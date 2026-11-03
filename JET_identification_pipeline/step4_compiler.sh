#!/bin/bash
set -euo pipefail

#### NeonDisco-TEAL — Step 4: Compile pVACbind results across sizes 8-11
#### Usage: step4_compiler.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR>
####
#### Combines filtered.tsv from pVACbind sizes 8-11 for one sample into:
####   <SAMPLE_OUTPUT_DIR>/<SAMPLE_NAME>_final_combined_data_TEAL.tsv
####
#### Keeps ALL pVACbind columns verbatim. Adds three columns:
####   Sample  — prepended at column 1
####   Size    — prepended at column 2
####   Junction — appended as LAST column (includes width from ids.txt)
####
#### Junction format per entry: chr:pos:strand:elem>chr:pos:strand:elem[width=N]
#### Multiple junctions for the same peptide separated by ;
####
#### Why width is included:
####   width=N tells how many amino acids span the exon-TE boundary.
####   It describes how much sequence was extracted around the breakpoint
####   and is relevant for reading frame interpretation and understanding
####   how much of the peptide comes from the TE vs the exon.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 2 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2

require_dir "${OUT_SAMPLE_DIR}" "Sample output directory"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"
RESULT_FILE="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}_final_combined_data_TEAL.tsv"
TMP_DIR_LOCAL=$(mktemp -d)
trap 'rm -rf "${TMP_DIR_LOCAL}"' EXIT

log_step "STEP 4 — Compiling pVACbind results: ${SAMPLE_NAME}"
log_info "Output file: ${RESULT_FILE}"

# =============================================================================
#  PRE-FLIGHT: CONFIRM AT LEAST ONE SIZE HAS BOTH REQUIRED FILES
# =============================================================================
FOUND_ANY=false
for SIZE in 8 9 10 11; do
    FILTERED="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.filtered.tsv"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"
    if [ -f "${FILTERED}" ] && [ -f "${IDSFILE}" ]; then
        FOUND_ANY=true
        break
    fi
done
[ "${FOUND_ANY}" = true ] || die "No pVACbind filtered.tsv + ids.txt found for ${SAMPLE_NAME}"

# =============================================================================
#  BUILD COMBINED OUTPUT
#
#  Column order:
#    1       : Sample          (added)
#    2       : Size            (added)
#    3..N    : all pVACbind filtered.tsv columns verbatim
#    last    : Junction        (added — junction coords + width from ids.txt)
#
#  Junction format:
#    chr19:44818637:+:ENST00000270233>chr19:44818713:+:(CCT)n[width=13]
#    Multiple junctions for the same peptide ID are ;-separated
#
#  ids.txt Name field format per combo (/-separated):
#    Order;ORFn;JunctionCoords;width=N
#    e.g. R1R2;ORF1;chr19:...>chr19:...;width=13
#    → we keep JunctionCoords[width=N] and deduplicate by JunctionCoords
# =============================================================================
HEADER_WRITTEN=false
> "${RESULT_FILE}"

for SIZE in 8 9 10 11; do
    FILTERED="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.filtered.tsv"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"

    if [ ! -f "${FILTERED}" ] || [ ! -f "${IDSFILE}" ]; then
        log_warn "  [size=${SIZE}] Missing filtered.tsv or ids.txt — skipping" >&2
        continue
    fi

    NROWS=$(tail -n +2 "${FILTERED}" | wc -l)
    if [ "${NROWS}" -eq 0 ]; then
        log_info "  [size=${SIZE}] 0 rows in filtered.tsv — skipping" >&2
        continue
    fi

    # Build Index → Junction[width] lookup
    # ids.txt col2 format: R1R2;ORF1;JunctionCoords;width=N / R1R2;ORF1;...
    # Extract JunctionCoords (part3) and width=N (part4) per unique junction
    awk -F'\t' '{
        n=split($2, combos, "/")
        delete seen; out=""
        for(i=1; i<=n; i++){
            split(combos[i], parts, ";")
            junc  = parts[3]          # junction coordinates
            width = parts[4]          # width=N  (may be empty for old files)
            key   = junc
            if(!(key in seen)){
                seen[key] = 1
                # Append [width=N] if present
                entry = (width != "") ? junc"["width"]" : junc
                out = (out == "") ? entry : out";"entry
            }
        }
        print $1"\t"out
    }' "${IDSFILE}" > "${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv"

    # Write header once — from first size file encountered
    # Column order: Sample | Size | <all pVACbind cols> | Junction
    if [ "${HEADER_WRITTEN}" = false ]; then
        PVAC_HEADER=$(head -1 "${FILTERED}")
        NCOLS=$(echo "${PVAC_HEADER}" | awk -F'\t' '{print NF}')
        echo -e "Sample\tSize\t${PVAC_HEADER}\tJunction" >> "${RESULT_FILE}"
        HEADER_WRITTEN=true
        log_info "  Header written: 2 prepended + ${NCOLS} pVACbind + 1 Junction = $((NCOLS+3)) total columns"
    fi

    # Append data rows
    awk -F'\t' \
        -v sample="${SAMPLE_NAME}" \
        -v sz="${SIZE}" \
        -v lookup="${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv" '
    BEGIN{
        while((getline line < lookup) > 0){
            split(line, a, "\t")
            junc[a[1]] = a[2]
        }
    }
    NR==1{ next }   # skip header already written
    {
        idx = $1
        j   = (idx in junc) ? junc[idx] : "NA"
        # Print: Sample | Size | all pVACbind cols | Junction
        printf "%s\t%s\t%s\t%s\n", sample, sz, $0, j
    }' "${FILTERED}" >> "${RESULT_FILE}"

    log_ok "  [size=${SIZE}] ${NROWS} rows written"
done

# =============================================================================
#  SUMMARY
# =============================================================================
TOTAL=$(tail -n +2 "${RESULT_FILE}" | wc -l)
TOTAL_COLS=$(head -1 "${RESULT_FILE}" | awk -F'\t' '{print NF}')

echo ""
log_info "── Summary for ${SAMPLE_NAME} ──────────────────────────"
log_info "  Total rows written     : ${TOTAL}"
log_info "  Total columns          : ${TOTAL_COLS}"
log_info "  Per-size breakdown:"
tail -n +2 "${RESULT_FILE}" | awk -F'\t' '{count[$2]++} END{
    for(s in count) print "    size"s": "count[s]
}' | sort
log_info "────────────────────────────────────────────────────────"
echo ""
log_ok "Step 4 complete. Results: ${RESULT_FILE}"
