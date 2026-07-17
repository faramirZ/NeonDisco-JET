#!/bin/bash
set -euo pipefail

#### JET Pipeline — Step 1: STAR alignment
#### Usage: step1_alignment.sh <SAMPLE_NAME> <R1.fastq.gz> <R2.fastq.gz>
#### Requires env vars (exported by master): STAR_BIN, SAMTOOLS_BIN,
####   GENOME_FASTA, GTF_FILE, STAR_INDEX_DIR, OUTPUTS_DIR, TMP_DIR,
####   READ_LENGTH, THREADS, LIB_DIR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${LIB_DIR:-$SCRIPT_DIR/lib}/common.sh"

[ "$#" -eq 3 ] || die "Usage: $(basename "$0") <SAMPLE_NAME> <R1.fastq.gz> <R2.fastq.gz>"
SAMPLE_NAME=$1
R1_GZ=$2
R2_GZ=$3

: "${STAR_BIN:?STAR_BIN not set}"
: "${SAMTOOLS_BIN:?SAMTOOLS_BIN not set}"
: "${GENOME_FASTA:?GENOME_FASTA not set}"
: "${GTF_FILE:?GTF_FILE not set}"
: "${STAR_INDEX_DIR:?STAR_INDEX_DIR not set}"
: "${OUTPUTS_DIR:?OUTPUTS_DIR not set}"
: "${TMP_DIR:?TMP_DIR not set}"
: "${READ_LENGTH:?READ_LENGTH not set}"
: "${THREADS:?THREADS not set}"

require_file "${R1_GZ}" "R1 fastq"
require_file "${R2_GZ}" "R2 fastq"
require_file "${STAR_BIN}/STAR" "STAR binary"

DAY=$(date +"%Y%m%d")
OUT_SAMPLE_DIR="${OUTPUTS_DIR}/${SAMPLE_NAME}_${DAY}"
mkdir -p "${OUT_SAMPLE_DIR}"
PREFIX="${OUT_SAMPLE_DIR}/${SAMPLE_NAME}"

R1=${R1_GZ%.gz}
R2=${R2_GZ%.gz}
[ -f "${R1}" ] || gunzip -k -c "${R1_GZ}" > "${R1}"
[ -f "${R2}" ] || gunzip -k -c "${R2_GZ}" > "${R2}"

SAM_FILE="${PREFIX}_Chimeric.out.sam"
BAM_FILE="${PREFIX}_Aligned.sortedByCoord.out.bam"
CHIM_BAM="${PREFIX}_Chimeric.out.bam"
CHIM_SORT_BAM="${PREFIX}_Chimeric.out.sort.bam"

log_step "STEP 1 — STAR alignment: ${SAMPLE_NAME}"
log_info "R1: ${R1}"
log_info "R2: ${R2}"
log_info "Output prefix: ${PREFIX}"

if [ ! -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
    log_info "STAR index not found — building genome index (this takes ~40 min)..."
    mkdir -p "${STAR_INDEX_DIR}"
    run_cmd "STAR genome indexing" \
        "${STAR_BIN}/STAR" \
            --runMode genomeGenerate \
            --genomeDir "${STAR_INDEX_DIR}" \
            --genomeFastaFiles "${GENOME_FASTA}" \
            --sjdbGTFfile "${GTF_FILE}" \
            --runThreadN "${THREADS}"
    log_ok "Genome index built successfully"
fi

run_cmd "STAR alignment (${SAMPLE_NAME})" \
    "${STAR_BIN}/STAR" \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX_DIR}" \
        --sjdbGTFfile "${GTF_FILE}" \
        --sjdbOverhang "${READ_LENGTH}" \
        --readFilesIn "${R1}" "${R2}" \
        --outFileNamePrefix "${PREFIX}_" \
        --outTmpDir "${TMP_DIR}/STAR_${SAMPLE_NAME}_${DAY}_$$" \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --bamRemoveDuplicatesType UniqueIdentical \
        --outFilterMismatchNoverLmax 0.04 \
        --outMultimapperOrder Random \
        --outFilterMultimapNmax 1000 \
        --winAnchorMultimapNmax 1000 \
        --chimOutType WithinBAM \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10

run_cmd "samtools view (chimeric)" \
    bash -c "'${SAMTOOLS_BIN}/samtools' view -@ '${THREADS}' -b '${SAM_FILE}' > '${CHIM_BAM}'"
run_cmd "samtools sort (chimeric)" \
    "${SAMTOOLS_BIN}/samtools" sort -@ "${THREADS}" -o "${CHIM_SORT_BAM}" -O bam "${CHIM_BAM}"
run_cmd "samtools index (chimeric)" "${SAMTOOLS_BIN}/samtools" index "${CHIM_SORT_BAM}"
run_cmd "samtools index (aligned)"  "${SAMTOOLS_BIN}/samtools" index "${BAM_FILE}"

log_ok "Step 1 complete for ${SAMPLE_NAME}. Output dir: ${OUT_SAMPLE_DIR}"
echo "${OUT_SAMPLE_DIR}"   # printed last on stdout so the caller can capture it
