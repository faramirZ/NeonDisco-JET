#!/bin/bash
set -euo pipefail

#### JET Pipeline — Step 4: Aggregate binders across peptide sizes 8-11
#### Usage: aggregate_binders.sh <SAMPLE_NAME> <SAMPLE_OUTPUT_DIR>

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

log_step "STEP 4 — Aggregating binders: ${SAMPLE_NAME}"

# Pre-check outside the subshell — FOUND_ANY set inside a pipe block
# runs in a subshell and cannot propagate back to the parent shell,
# so we check file existence here before entering the pipe.
FOUND_ANY=false
for SIZE in 8 9 10 11; do
    NETFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.netmhcpan4.txt"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"
    if [ -f "${NETFILE}" ] && [ -f "${IDSFILE}" ]; then
        FOUND_ANY=true
        break
    fi
done
[ "${FOUND_ANY}" = true ] || die "No netMHCpan output files found for ${SAMPLE_NAME}"

# Build the combined binder table across sizes 8-11
{
echo -e "Size\tID\tPeptide\tMinRank\tNB\tJunction"

for SIZE in 8 9 10 11; do
    NETFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.netmhcpan4.txt"
    IDSFILE="${PREFIX}_Fusions_chim2.junc2.size${SIZE}.ids.txt"
    [ -f "${NETFILE}" ] && [ -f "${IDSFILE}" ] || { log_warn "Skipping size ${SIZE}: missing netMHCpan or ids file"; continue; }

    awk -F'\t' '{
        n=split($2, combos, "/")
        delete seen
        out=""
        for(i=1;i<=n;i++){
            split(combos[i], parts, ";")
            junc = parts[3]
            if(!(junc in seen)){ seen[junc]=1; out = (out=="") ? junc : out";"junc }
        }
        print $1"\t"out
    }' "${IDSFILE}" > "${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv"

    # Detect number of alleles from netMHCpan header line (line 3)
    # Column layout: Pos(1) Peptide(2) ID(3) [core icore EL_score EL_rank]*N Ave NB
    # Total cols = 4*N + 5  →  N = (NF - 5) / 4
    # Rank column for allele i = 3 + 4*i  (e.g. allele1=col7, allele2=col11, ...)
    N_ALLELES=$(awk -F'\t' 'NR==3{ print int((NF-5)/4); exit }' "${NETFILE}")

    # Log allele count and rank column positions for this size file
    log_info "  [size=${SIZE}] Detected ${N_ALLELES} allele(s) in netMHCpan output"
    RANK_COL_INFO=""
    for (( i=1; i<=N_ALLELES; i++ )); do
        COL=$(( 3 + 4*i ))
        RANK_COL_INFO="${RANK_COL_INFO} allele${i}→col${COL}"
    done
    log_info "  [size=${SIZE}] Rank columns:${RANK_COL_INFO}"

    awk -F'\t' -v sz="${SIZE}" -v n_alleles="${N_ALLELES}" -v lookup="${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv" '
    BEGIN{
        while((getline line < lookup) > 0){ split(line,a,"\t"); junc[a[1]]=a[2] }
        for (i=1; i<=n_alleles; i++) rank_col[i] = 3 + 4*i
    }
    NR>3 {
        minrank = $(rank_col[1])
        for (i=2; i<=n_alleles; i++)
            if ($(rank_col[i]) < minrank) minrank = $(rank_col[i])
        if (minrank < 2) {
            j = ($3 in junc) ? junc[$3] : "NA"
            printf "%s\t%s\t%s\t%.4f\t%s\t%s\n", sz, $3, $2, minrank, $NF, j
        }
    }' "${NETFILE}"
done
} | awk -F'\t' 'NR==1{print;next} !seen[$2"\t"$3]++' > "${RESULT_FILE}"

TOTAL=$(tail -n +2 "${RESULT_FILE}" | wc -l)
STRONG=$(tail -n +2 "${RESULT_FILE}" | awk -F'\t' '$4<0.5' | wc -l)
WEAK=$((TOTAL - STRONG))
UNIQ_JUNCTIONS=$(tail -n +2 "${RESULT_FILE}" | awk -F'\t' '{print $6}' | tr ';' '\n' | sort -u | wc -l)

log_info "No of Binders (rank<2, sizes 8-11): ${TOTAL}"
log_info "No of Strong Binders (rank<0.5)  : ${STRONG}"
log_info "No of Weak Binders (0.5<=rank<2) : ${WEAK}"
log_info "No of unique JET junctions       : ${UNIQ_JUNCTIONS}"
log_info "Per-size breakdown:"
tail -n +2 "${RESULT_FILE}" | awk -F'\t' '{cat=($4<0.5)?"Strong":"Weak"; key=$1"\t"cat; count[key]++} END{for(k in count) print "          "k": "count[k]}' | sort

log_ok "Step 4 complete. Results: ${RESULT_FILE}"
