#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "Usage: $0 R1 R2 SAMPLE THREADS STAR_INDEX GTF SJDB_OVERHANG" >&2
    exit 1
fi

r1=$1
r2=$2
sample=$3
threads=$4
star_index=$5
gtf=$6
sjdb_overhang=$7

for input in "$r1" "$r2" "$gtf"; do
    [[ -r "$input" ]] || {
        echo "Cannot read input: $input" >&2
        exit 1
    }
done

read_args=()
if [[ "$r1" == *.gz && "$r2" == *.gz ]]; then
    read_args=(--readFilesCommand zcat)
elif [[ "$r1" == *.gz || "$r2" == *.gz ]]; then
    echo "Both FASTQ files must use the same compression format." >&2
    exit 1
fi

# Decompress a gzip GTF into this task's working directory.
if [[ "$gtf" == *.gz ]]; then
    plain_gtf=$(mktemp ./star-annotation.XXXXXX.gtf)
    trap 'rm -f -- "$plain_gtf"' EXIT
    gzip -dc -- "$gtf" > "$plain_gtf"
    gtf=$plain_gtf
fi

STAR \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --runThreadN "$threads" \
    --genomeDir "$star_index" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang "$sjdb_overhang" \
    --readFilesIn "$r1" "$r2" \
    "${read_args[@]}" \
    --outFileNamePrefix "${sample}_" \
    --outTmpDir "${PWD}/STAR_tmp_${sample}_$$" \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMismatchNoverLmax 0.04 \
    --outMultimapperOrder Random \
    --outFilterMultimapNmax 1000 \
    --winAnchorMultimapNmax 1000 \
    --chimOutType Junctions WithinBAM \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 10
    