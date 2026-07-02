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

    awk -F'\t' -v sz="${SIZE}" -v lookup="${TMP_DIR_LOCAL}/lookup_size${SIZE}.tsv" '
    BEGIN{ while((getline line < lookup) > 0){ split(line,a,"\t"); junc[a[1]]=a[2] } }
    NR>3 {
        minrank=$7
        if($11<minrank) minrank=$11
        if($15<minrank) minrank=$15
        if($19<minrank) minrank=$19
        if($23<minrank) minrank=$23
        if($27<minrank) minrank=$27
        if(minrank<2){
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
