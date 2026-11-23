(
base=/home/faramir/data/neondisco-jet/outputData
old="$base/379T_20261104"
new="$base/379T_nextflow_steps1to/pVACbind/379T"

if [[ ! -d "$new" ]]; then
    new="$base/379T_nextflow_steps1to3/pVACbind/379T"
fi
for dir in "$old" "$new"; do
    [[ -d "$dir" ]] || { echo "Missing directory: $dir" >&2; exit 2; }
done
echo "Bash:     $old"
echo "Nextflow: $new"

matched=0
for size in 8 9 10 11; do
    for report in all_epitopes filtered all_epitopes.aggregated; do
        file="pvacbind_size$size/MHC_Class_I/379T.MHC_I.$report.tsv"

        if [[ ! -f "$old/$file" || ! -f "$new/$file" ]]; then
            [[ -f "$old/$file" ]] || echo "MISSING BASH      $file"
            [[ -f "$new/$file" ]] || echo "MISSING NEXTFLOW  $file"
        elif cmp -s -- "$old/$file" "$new/$file"; then
            echo "IDENTICAL         $file"
            ((matched+=1))
        elif [[ $? -eq 1 ]]; then
            echo "DIFFERENT         $file"
        else
            echo "READ ERROR        $file"
        fi
    done
done
echo "$matched of 12 reports are byte-for-byte identical."
)
