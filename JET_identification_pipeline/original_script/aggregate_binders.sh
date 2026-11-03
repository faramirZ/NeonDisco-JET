#!/bin/bash

outdir="${SAMPLE_OUTDIR:-/home/faramir/data/neondisco-jet/outputData/379T_20260615}"
prefix="${SAMPLE_PREFIX:-379T}"
result_file="${outdir}/all_sizes_binders.tsv"
tmp_dir=$(mktemp -d)

{
echo -e "Size\tID\tPeptide\tMinRank\tNB\tJunction"

for size in 8 9 10 11; do
  netfile="${outdir}/${prefix}_Fusions_chim2.junc2.size${size}.netmhcpan4.txt"
  idsfile="${outdir}/${prefix}_Fusions_chim2.junc2.size${size}.ids.txt"
  [ -f "$netfile" ] || continue
  [ -f "$idsfile" ] || continue

  # Build ID -> Junction lookup for this size.
  # ids.txt Name format: Order;ORFn;JunctionId;width=N  (multiple combos joined by "/")
  # We extract only the JunctionId (3rd ";" field) for each combo, dedupe, join by ";"
  awk -F'\t' '{
    n=split($2, combos, "/")
    delete seen
    out=""
    for(i=1;i<=n;i++){
      split(combos[i], parts, ";")
      junc = parts[3]
      if(!(junc in seen)){
        seen[junc]=1
        out = (out=="") ? junc : out";"junc
      }
    }
    print $1"\t"out
  }' "$idsfile" > "${tmp_dir}/lookup_size${size}.tsv"

  awk -F'\t' -v sz="$size" -v lookup="${tmp_dir}/lookup_size${size}.tsv" '
  BEGIN{
    while((getline line < lookup) > 0){
      split(line, a, "\t")
      junc[a[1]] = a[2]
    }
  }
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
  }' "$netfile"
done
} | awk -F'\t' 'NR==1{print;next} !seen[$2"\t"$3]++' > "${result_file}"

rm -rf "${tmp_dir}"

total=$(tail -n +2 "${result_file}" | wc -l)
strong=$(tail -n +2 "${result_file}" | awk -F'\t' '$4<0.5' | wc -l)
weak=$((total - strong))
uniq_junctions=$(tail -n +2 "${result_file}" | awk -F'\t' '{print $6}' | tr ';' '\n' | sort -u | wc -l)

echo "No of Binders (rank<2, sizes 8-11): ${total}"
echo "No of Strong Binders (rank<0.5): ${strong}"
echo "No of Weak Binders (0.5<=rank<2): ${weak}"
echo "No of unique JET junctions represented: ${uniq_junctions}"
echo ""
echo "Per-size breakdown:"
tail -n +2 "${result_file}" | awk -F'\t' '{cat=($4<0.5)?"Strong":"Weak"; key=$1"\t"cat; count[key]++} END{for(k in count) print k": "count[k]}' | sort
echo ""
echo "Full binder list saved to: ${result_file}"
echo "(columns: Size, ID, Peptide, MinRank, NB, Junction)"
