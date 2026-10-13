#!/bin/bash
set -euo pipefail

#### NeonDisco-TEAL — Step 4: Aggregate pVACbind binders across peptide sizes 8-11
#### Usage: aggregate_binders.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR>
####
#### Reads pVACbind all_epitopes.tsv per size, maps peptides back to JET
#### junctions via per-size ids.txt, filters by IC50 threshold, and writes
#### a single combined all_sizes_binders.tsv for the sample.
####
#### Binder thresholds (IC50 nM):
####   Strong binder : Best IC50 <  50 nM
####   Weak binder   : Best IC50 < 500 nM
####
#### Output columns:
####   Size | Index | Peptide | HLA_Allele | Best_IC50 | Median_IC50 |
####   Best_Percentile | Best_Method | Junction

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 2 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR>"
SAMPLE_NAME=$1
OUT_SAMPLE_DIR=$2

require_dir "${OUT_SAMPLE_DIR}" "Sample output directory"

PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"
RESULT_FILE="${OUT_SAMPLE_DIR}/all_sizes_binders.tsv"
TMP_DIR_LOCAL=$(mktemp -d)
trap 'rm -rf "${TMP_DIR_LOCAL}"' EXIT

log_step "STEP 4 — Aggregating pVACbind binders: ${SAMPLE_NAME}"

# =============================================================================
#  PRE-FLIGHT: CONFIRM AT LEAST ONE SIZE HAS BOTH REQUIRED FILES
# =============================================================================
FOUND_ANY=false
for SIZE in 8 9 10 11; do
    PVAC_FILE="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.all_epitopes.tsv"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"
    if [ -f "${PVAC_FILE}" ] && [ -f "${IDSFILE}" ]; then
        FOUND_ANY=true
        break
    fi
done
[ "${FOUND_ANY}" = true ] || die "No pVACbind output files found for ${SAMPLE_NAME}"

# =============================================================================
#  BUILD COMBINED BINDER TABLE
#
#  all_epitopes.tsv column layout (confirmed from test run):
#    1  Index                    — FASTA entry ID (matches ids.txt col 1)
#    2  HLA Allele
#    3  Sub-peptide Position
#    4  Epitope Seq              — peptide amino acid sequence
#    5  Median IC50 Score        — median across all algorithms
#    6  Best IC50 Score          — best (lowest) IC50 across all algorithms
#    7  Best IC50 Score Method   — which algorithm gave the best IC50
#    8  Median Percentile
#    9  Best Percentile          — best percentile across all algorithms
#   10  Best Percentile Method
#   ... per-algorithm columns follow
#
#  Each row is one peptide × one HLA allele combination.
#  We dedup within each size by Index, keeping the allele row with lowest
#  Best IC50, then filter to IC50 < 500 nM (binders).
# =============================================================================
{
echo -e "Size\tIndex\tPeptide\tHLA_Allele\tBest_IC50\tMedian_IC50\tBest_Percentile\tBest_Method\tJunction"

for SIZE in 8 9 10 11; do
    PVAC_FILE="${OUT_SAMPLE_DIR}/pvacbind_size${SIZE}/MHC_Class_I/${SAMPLE_NAME}.MHC_I.all_epitopes.tsv"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"

    [ -f "${PVAC_FILE}" ] && [ -f "${IDSFILE}" ] \
        || { log_warn "Skipping size ${SIZE}: missing pVACbind output or ids file" >&2; continue; }

    # Build Index → Junction lookup from ids.txt
    # ids.txt Name field format: Order;ORFn;JunctionId;width=N (multiple combos joined by /)
    awk -F'\t' '{
        n=split($2, combos, "/")
        delete seen
        out=""
        for(i=1; i<=n; i++){
            split(combos[i], parts, ";")
            junc = parts[3]
            if(!(junc in seen)){ seen[junc]=1; out = (out=="") ? junc : out";"junc }
        }
        print $1"\t"out
    }' "${IDSFILE}" > "${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv"

    NLINES=$(tail -n +2 "${PVAC_FILE}" | wc -l)
    log_info "  [size=${SIZE}] ${NLINES} epitope-allele rows in pVACbind output" >&2

    # Parse all_epitopes.tsv, dedup by Index (keep best IC50 across alleles),
    # filter to IC50 < 500 nM, map to junction, emit one row per passing peptide
    awk -F'\t' -v sz="${SIZE}" -v lookup="${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv" '
    BEGIN{
        while((getline line < lookup) > 0){
            split(line, a, "\t")
            junc[a[1]] = a[2]
        }
    }
    NR==1{ next }   # skip header
    {
        idx     = $1
        allele  = $2
        peptide = $4
        median  = $5
        best    = $6
        method  = $7
        pct     = $9

        # Skip rows with missing or non-numeric IC50
        if(best == "NA" || best == "" || best+0 != best+0) next

        # Only consider binders (IC50 < 500 nM)
        if(best+0 >= 500) next

        # Keep row with lowest Best IC50 per Index (best allele wins)
        if(!(idx in min_ic50) || best+0 < min_ic50[idx]+0){
            min_ic50[idx]    = best
            min_median[idx]  = median
            min_allele[idx]  = allele
            min_peptide[idx] = peptide
            min_pct[idx]     = pct
            min_method[idx]  = method
        }
    }
    END{
        for(idx in min_ic50){
            j = (idx in junc) ? junc[idx] : "NA"
            printf "%s\t%s\t%s\t%s\t%.4f\t%.4f\t%s\t%s\t%s\n",
                sz,
                idx,
                min_peptide[idx],
                min_allele[idx],
                min_ic50[idx]+0,
                min_median[idx]+0,
                min_pct[idx],
                min_method[idx],
                j
        }
    }' "${PVAC_FILE}"

done

# Outer dedup: remove identical Size+Index+Peptide rows that slipped through
} | awk -F'\t' 'NR==1{print;next} !seen[$1"\t"$2"\t"$3]++' > "${RESULT_FILE}"

# =============================================================================
#  SUMMARY COUNTS
#  Strong binder : Best_IC50 (col 5) <  50 nM
#  Binder        : Best_IC50 (col 5) < 500 nM  (all rows qualify — already filtered above)
#  Unique junctions from col 9
# =============================================================================
TOTAL=$(tail -n +2 "${RESULT_FILE}" | wc -l)
STRONG=$(tail -n +2 "${RESULT_FILE}" | awk -F'\t' '$5<50'  | wc -l)
WEAK=$((TOTAL - STRONG))
UNIQ_JUNCTIONS=$(tail -n +2 "${RESULT_FILE}" | awk -F'\t' '{print $9}' | tr ';' '\n' | sort -u | wc -l)

log_info "No of Binders (IC50<500nM, sizes 8-11) : ${TOTAL}"
log_info "No of Strong Binders (IC50<50nM)       : ${STRONG}"
log_info "No of Weak Binders (50<=IC50<500nM)    : ${WEAK}"
log_info "No of unique JET junctions             : ${UNIQ_JUNCTIONS}"
log_info "Per-size breakdown:"
tail -n +2 "${RESULT_FILE}" \
    | awk -F'\t' '{cat=($5<50)?"Strong":"Weak"; key=$1"\t"cat; count[key]++} END{for(k in count) print "          "k": "count[k]}' \
    | sort

log_ok "Step 4 complete. Results: ${RESULT_FILE}"
