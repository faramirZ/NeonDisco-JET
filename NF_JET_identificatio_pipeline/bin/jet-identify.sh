#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 8 ]]; then
    echo "Usage: $0 SAMPLE CHIMERIC JUNCTIONS STAR_LOG REPEATS GTF R_SCRIPT GENOME" >&2
    exit 1
fi

sample=$1
chimeric=$2
junctions=$3
star_log=$4
repeats=$5
gtf=$6
r_script=$7
genome=$8

rscript_bin="${RSCRIPT_BIN:-Rscript}"

command -v "$rscript_bin" >/dev/null || {
    echo "Rscript executable not found: $rscript_bin" >&2
    exit 1
}

for input in "$chimeric" "$junctions" "$star_log" "$repeats" "$gtf" "$r_script"; do
    [[ -r "$input" ]] || {
        echo "Cannot read input: $input" >&2
        exit 1
    }
done

libsize=$(awk -F '|' '
    $1 ~ /^[[:space:]]*Uniquely mapped reads number[[:space:]]*$/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
    }
' "$star_log")

if [[ ! "$libsize" =~ ^[0-9]+$ || ! "$libsize" =~ [1-9] ]]; then
    echo "Expected a positive uniquely mapped read count in $star_log; found: $libsize" >&2
    exit 1
fi

echo "Sample: $sample; uniquely mapped reads: $libsize"

for size in 8 9 10 11; do
    "$rscript_bin" --vanilla "$r_script" \
        --chimeric "$chimeric" \
        --junction "$junctions" \
        --genome "$genome" \
        --size "$size" \
        --libsize "$libsize" \
        --repeats "$repeats" \
        --gtf "$gtf" \
        --prefix "$sample" \
        --verbose
done
