#!/bin/bash

#### JET Pipeline — Output Integrity Check
#### Compares a new pipeline run against a reference run using md5sum

OLD_DIR=/home/faramir/data/neondisco-jet/outputData/379T_20260615
NEW_DIR=$(ls -td /home/faramir/data/neondisco-jet/outputData/379T_2* 2>/dev/null | grep -v "20260615" | head -1)
SAMPLE=379T
REPORT_DIR=/home/faramir/data/neondisco-jet/outputData/integrity_check
mkdir -p "${REPORT_DIR}"
REPORT="${REPORT_DIR}/integrity_report_$(date +%Y%m%d_%H%M%S).txt"

echo "========================================" | tee "${REPORT}"
echo "JET Pipeline — Integrity Check Report"   | tee -a "${REPORT}"
echo "Reference : ${OLD_DIR}"                  | tee -a "${REPORT}"
echo "New run   : ${NEW_DIR}"                  | tee -a "${REPORT}"
echo "Generated : $(date)"                     | tee -a "${REPORT}"
echo "========================================" | tee -a "${REPORT}"

if [ -z "${NEW_DIR}" ]; then
    echo "ERROR: Could not find a new run directory. Exiting." | tee -a "${REPORT}"
    exit 1
fi

#-------------------------------------------------------
# Define file groups to check
#-------------------------------------------------------

declare -A FILE_GROUPS
FILE_GROUPS["STAR_outputs"]="${SAMPLE}_Log.final.out ${SAMPLE}_SJ.out.tab ${SAMPLE}_ReadsPerGene.out.tab"
FILE_GROUPS["BAM_indexes"]="${SAMPLE}_Aligned.sortedByCoord.out.bam.bai ${SAMPLE}_Chimeric.out.sort.bam.bai"
FILE_GROUPS["JET_annotation"]="${SAMPLE}_Chimeric.out.annotatedJET.txt"
FILE_GROUPS["size8"]="${SAMPLE}_Fusions_chim2.junc2.size8.fasta ${SAMPLE}_Fusions_chim2.junc2.size8.genomic.txt ${SAMPLE}_Fusions_chim2.junc2.size8.ids.txt ${SAMPLE}_Fusions_chim2.junc2.size8.netmhcpan4.txt"
FILE_GROUPS["size9"]="${SAMPLE}_Fusions_chim2.junc2.size9.fasta ${SAMPLE}_Fusions_chim2.junc2.size9.genomic.txt ${SAMPLE}_Fusions_chim2.junc2.size9.ids.txt ${SAMPLE}_Fusions_chim2.junc2.size9.netmhcpan4.txt"
FILE_GROUPS["size10"]="${SAMPLE}_Fusions_chim2.junc2.size10.fasta ${SAMPLE}_Fusions_chim2.junc2.size10.genomic.txt ${SAMPLE}_Fusions_chim2.junc2.size10.ids.txt ${SAMPLE}_Fusions_chim2.junc2.size10.netmhcpan4.txt"
FILE_GROUPS["size11"]="${SAMPLE}_Fusions_chim2.junc2.size11.fasta ${SAMPLE}_Fusions_chim2.junc2.size11.genomic.txt ${SAMPLE}_Fusions_chim2.junc2.size11.ids.txt ${SAMPLE}_Fusions_chim2.junc2.size11.netmhcpan4.txt"
FILE_GROUPS["aggregation"]="all_sizes_binders.tsv"

#-------------------------------------------------------
# Compare each file group
#-------------------------------------------------------

PASS=0
FAIL=0
SKIP=0
WARN=0

check_file() {
    local group="$1"
    local fname="$2"
    local old_path="${OLD_DIR}/${fname}"
    local new_path="${NEW_DIR}/${fname}"

    printf "  %-60s" "${fname}" | tee -a "${REPORT}"

    # Check existence
    if [ ! -f "${old_path}" ] && [ ! -f "${new_path}" ]; then
        printf "[SKIP] both missing\n" | tee -a "${REPORT}"
        SKIP=$((SKIP+1)); return
    fi
    if [ ! -f "${old_path}" ]; then
        printf "[WARN] only in new run\n" | tee -a "${REPORT}"
        WARN=$((WARN+1)); return
    fi
    if [ ! -f "${new_path}" ]; then
        printf "[FAIL] missing in new run\n" | tee -a "${REPORT}"
        FAIL=$((FAIL+1)); return
    fi

    # Line counts (more meaningful than md5 for text files
    # where timestamps or ordering may differ slightly)
    old_lines=$(wc -l < "${old_path}")
    new_lines=$(wc -l < "${new_path}")

    # md5sum
    old_md5=$(md5sum "${old_path}" | cut -d' ' -f1)
    new_md5=$(md5sum "${new_path}" | cut -d' ' -f1)

    if [ "${old_md5}" = "${new_md5}" ]; then
        printf "[PASS] identical (md5: %s, lines: %s)\n" "${old_md5:0:12}..." "${old_lines}" | tee -a "${REPORT}"
        PASS=$((PASS+1))
    elif [ "${old_lines}" = "${new_lines}" ]; then
        printf "[WARN] md5 differs but same line count (%s lines) — check content\n" "${old_lines}" | tee -a "${REPORT}"
        WARN=$((WARN+1))
    else
        printf "[FAIL] md5 differs AND line count differs (old: %s, new: %s lines)\n" "${old_lines}" "${new_lines}" | tee -a "${REPORT}"
        FAIL=$((FAIL+1))
    fi
}

for group in "STAR_outputs" "BAM_indexes" "JET_annotation" "size8" "size9" "size10" "size11" "aggregation"; do
    echo "" | tee -a "${REPORT}"
    echo "--- ${group} ---" | tee -a "${REPORT}"
    for fname in ${FILE_GROUPS[$group]}; do
        check_file "${group}" "${fname}"
    done
done

#-------------------------------------------------------
# Deep content check for key files
#-------------------------------------------------------

echo "" | tee -a "${REPORT}"
echo "--- Deep content check (key files) ---" | tee -a "${REPORT}"

# Compare binder counts from aggregation
for run_label in "OLD:${OLD_DIR}" "NEW:${NEW_DIR}"; do
    label="${run_label%%:*}"
    dir="${run_label##*:}"
    f="${dir}/all_sizes_binders.tsv"
    if [ -f "$f" ]; then
        total=$(tail -n +2 "$f" | wc -l)
        strong=$(tail -n +2 "$f" | awk -F'\t' '$4<0.5' | wc -l)
        weak=$((total-strong))
        uniq=$(tail -n +2 "$f" | awk -F'\t' '{print $6}' | tr ';' '\n' | sort -u | wc -l)
        echo "  ${label} binders: total=${total}  strong=${strong}  weak=${weak}  unique_junctions=${uniq}" | tee -a "${REPORT}"
    else
        echo "  ${label} all_sizes_binders.tsv not found" | tee -a "${REPORT}"
    fi
done

# Compare annotatedJET line counts
echo "" | tee -a "${REPORT}"
for run_label in "OLD:${OLD_DIR}" "NEW:${NEW_DIR}"; do
    label="${run_label%%:*}"
    dir="${run_label##*:}"
    f="${dir}/${SAMPLE}_Chimeric.out.annotatedJET.txt"
    if [ -f "$f" ]; then
        total=$(tail -n +2 "$f" | wc -l)
        jets=$(tail -n +2 "$f" | awk -F'\t' '$9=="TRUE"' | wc -l)
        echo "  ${label} annotatedJET: total_junctions=${total}  JETs_passing_filter=${jets}" | tee -a "${REPORT}"
    fi
done

# BAM file sizes (md5 will always differ due to timestamps in BAM header)
echo "" | tee -a "${REPORT}"
echo "  BAM file sizes (md5 not meaningful for BAMs — comparing size instead):" | tee -a "${REPORT}"
for fname in "${SAMPLE}_Aligned.sortedByCoord.out.bam" "${SAMPLE}_Chimeric.out.sort.bam"; do
    old_size=$(du -sh "${OLD_DIR}/${fname}" 2>/dev/null | cut -f1)
    new_size=$(du -sh "${NEW_DIR}/${fname}" 2>/dev/null | cut -f1)
    printf "  %-55s OLD: %-8s NEW: %s\n" "${fname}" "${old_size:-missing}" "${new_size:-missing}" | tee -a "${REPORT}"
done

#-------------------------------------------------------
# Summary
#-------------------------------------------------------

echo "" | tee -a "${REPORT}"
echo "========================================" | tee -a "${REPORT}"
echo "SUMMARY"                                  | tee -a "${REPORT}"
echo "  PASS : ${PASS}"                         | tee -a "${REPORT}"
echo "  WARN : ${WARN}  (check manually)"       | tee -a "${REPORT}"
echo "  FAIL : ${FAIL}"                         | tee -a "${REPORT}"
echo "  SKIP : ${SKIP}  (missing in both)"      | tee -a "${REPORT}"
echo "========================================" | tee -a "${REPORT}"
echo "Full report saved to: ${REPORT}"          | tee -a "${REPORT}"

[ "${FAIL}" -eq 0 ] && exit 0 || exit 1
