#!/usr/bin/env bash

READ1=$1
READ2=$2
SAMPLE_ID=$3
CORES=$4
INDEX=$5

# check if the input reads are gzipped or not, important for STAR command
case $READ1 in
    *.gz)
    # run STAR with zcat flag
        echo "Running STAR with zcat flag"
        STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_ID}-STAR_" \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMattributes All \
        --outFilterType BySJout \
        --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000
        ;;
    *.fastq | *.fq)
    # run STAR without zcat flag
        echo "Running STAR without zcat flag"
        STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --outFileNamePrefix "${SAMPLE_ID}-STAR_" \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMattributes All \
        --outFilterType BySJout \
        --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000
        ;;
    *)
        echo "Input reads are not in fastq or fastq.gz format"
        exit 1
        ;;
esac

