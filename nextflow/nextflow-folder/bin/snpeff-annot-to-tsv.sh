#!/usr/bin/env bash
# snpeff-annot-to-tsv.sh
# Extract and flatten SnpEff-annotated variants to primary TSV outputs.

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 allSingleSubsEditing.annot.vcf sample_name" >&2
    exit 1
fi

INPUT="$1"
SAMPLE_NAME="$2"

# Define the explicit static column mapping layout for SnpSift extractFields
# Field positions inside awk:
# $1=CHROM, $2=POS, $3=REF, $4=ALT, $5=ANN, $6=STRAND, $7=FREQ, $8=COV, $9=REPTYPE, $10=REPNAME, $11=REDIPORTAL, $12=EDITSTATUS
FIELDS="CHROM POS REF ALT ANN STRAND FREQ COV REPTYPE REPNAME REDIPORTAL EDITSTATUS"

echo "Extracting annotations to standardized streams..."

SnpSift extractFields "$INPUT" $FIELDS | awk -F'\t' -v prefix="${SAMPLE_NAME}" '
NR==1 {
    print $0 > (prefix "_single-subs-all-sites.tsv")
    print $0 > (prefix "_AG-subs-all-sites.tsv")
    print $0 > (prefix "_AG-subs-rediportal-known-sites.tsv")
}
NR>1 {
    # 1. Master Output (Already filtered for single-subs in VCF conversion)
    print $0 > (prefix "_single-subs-all-sites.tsv")

    # Define Canonical (A->G or T->C)
    is_canonical = (($3=="A" && $4=="G") || ($3=="T" && $4=="C"))

    # 2. Canonical Sites (all)
    if (is_canonical) {
        print $0 > (prefix "_AG-subs-all-sites.tsv")
    }

    # 3. Canonical Sites (REDIPortal known)
    if (is_canonical && $11 == "ed") {
        print $0 > (prefix "_AG-subs-rediportal-known-sites.tsv")
    }
}'

echo "Extraction complete. All primary TSV files generated."

# Ensure files exist and have content before proceeding to checks
for f in "${SAMPLE_NAME}"_*.tsv; do
    if [[ ! -s "$f" ]]; then
        echo "WARNING: File $f is empty. Skipping integrity checks for this sample." >&2
        exit 0 
    fi
done

echo "Performing automated integrity checks..."

# 1. Hierarchical Integrity Check
# We use local variables to ensure logic holds
MASTER_COUNT=$(grep -v "CHROM" "${SAMPLE_NAME}_single-subs-all-sites.tsv" | wc -l)
AG_ALL_COUNT=$(grep -v "CHROM" "${SAMPLE_NAME}_AG-subs-all-sites.tsv" | wc -l)
AG_KNOWN_COUNT=$(grep -v "CHROM" "${SAMPLE_NAME}_AG-subs-rediportal-known-sites.tsv" | wc -l)

if [ "$AG_ALL_COUNT" -gt "$MASTER_COUNT" ] || [ "$AG_KNOWN_COUNT" -gt "$AG_ALL_COUNT" ]; then
    echo "ERROR: Data hierarchy violation (Subset > Supersets)." >&2
    exit 1
fi

# 2. Canonical Leakage Check
# If this awk command finds anything, it exits with error code 1
awk -F'\t' 'NR>1 {
    if (!(($3=="A" && $4=="G") || ($3=="T" && $4=="C"))) {
        print "LEAKAGE: Non-canonical site found in AG file: " $1 ":" $2 " (" $3 ">" $4 ")"
        exit 1
    }
}' "${SAMPLE_NAME}_AG-subs-all-sites.tsv"

echo "Integrity checks passed. Pipeline finalized."

