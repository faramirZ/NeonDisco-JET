#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 SAMPLE FASTA HLA_ALLELES SIZE THREADS OUTDIR" >&2
    exit 1
fi

sample=$1
fasta=$2
alleles=$3
size=$4
threads=$5
outdir=$6

[[ -f "$fasta" && -r "$fasta" ]] || { echo "Cannot read $fasta" >&2; exit 1; }
[[ "$alleles" =~ HLA-[ABC] ]] || { echo "Invalid HLA list: $alleles" >&2; exit 1; }
alleles=$(sed -E 's/HLA-([ABC])([0-9])/HLA-\1*\2/g' <<< "$alleles")

mkdir -p "$outdir"
if [[ ! -s "$fasta" ]]; then
    echo "No peptides for $sample size $size" >&2
    exit 0
fi

extra=()
if [[ ${KEEP_TMP:-false} == true ]]; then
    extra+=(-k)
fi

pvacbind run "$fasta" "$sample" "$alleles" \
    NetMHCpan NetMHCpanEL MHCflurryEL BigMHC_IM PRIME DeepImmuno \
    "$outdir" \
    --iedb-install-directory /opt/iedb \
    --allele-specific-binding-thresholds \
    --netmhc-stab \
    -e1 "$size" -t "$threads" "${extra[@]}"
