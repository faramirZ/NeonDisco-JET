#!/usr/bin/env bash
# snpeff-annot-to-tsv.sh
# Extract SnpEff-annotated VCF to TSV format using SnpSift
# Produces two outputs:
#   1. All known editing sites (FREQ > 0.1, REDIPortal confirmed)
#   2. Known genic sites only (excludes intergenic)
#
# Usage: ./snpeff-annot-to-tsv.sh input.ann.vcf output_prefix

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.ann.vcf output_prefix" >&2
    exit 1
fi

INPUT="$1"
PREFIX="$2"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input file '$INPUT' not found" >&2
    exit 1
fi

FIELDS="CHROM POS REF ALT \
    ANN[0].GENE \
    ANN[0].FEATUREID \
    ANN[0].EFFECT \
    ANN[0].IMPACT \
    ANN[0].BIOTYPE \
    ANN[0].HGVS_C \
    STRAND \
    FREQ \
    COV \
    REPTYPE \
    REPNAME \
    EDITSTATUS"

# 1. All known sites with FREQ > 0.1
echo "Extracting all known sites (FREQ > 0.1)..."
SnpSift filter "FREQ > 0.1 & REDIPORTAL == 'ed'" "$INPUT" | \
SnpSift extractFields - $FIELDS \
    > "${PREFIX}.known_sites_freq10pct.tsv"

echo "  -> ${PREFIX}.known_sites_freq10pct.tsv ($(grep -c $'^[^#]' "${PREFIX}.known_sites_freq10pct.tsv") sites)"

# 2. Known genic sites only (intronic, UTR, exonic — excludes intergenic)
echo "Extracting known genic sites (FREQ > 0.1)..."
SnpSift filter "FREQ > 0.1 & REDIPORTAL == 'ed' & !(ANN[0].EFFECT has 'intergenic_region')" "$INPUT" | \
SnpSift extractFields - $FIELDS \
    > "${PREFIX}.known_sites_freq10pct_genic.tsv"

echo "  -> ${PREFIX}.known_sites_freq10pct_genic.tsv ($(grep -c $'^[^#]' "${PREFIX}.known_sites_freq10pct_genic.tsv") sites)"


# 3. Collapse gene-level annotations (one line per gene, with all effects concatenated)
echo "Collapsing gene-level annotations..."

# check if the input file exists and is not empty
if [[ ! -s "${PREFIX}.known_sites_freq10pct_genic.tsv" ]]; then
    echo "Error: input file '${PREFIX}.known_sites_freq10pct_genic.tsv' not found or empty" >&2
    exit 1
fi

(awk -F'\t' 'NR>1 && $5 != "." {
    gene=$5
    freq=$12
    effect=$7
    biotype=$9
    count[gene]++
    if (!minfreq[gene] || freq < minfreq[gene]) minfreq[gene]=freq
    if (freq > maxfreq[gene]) maxfreq[gene]=freq
    if (!seen_effect[gene,effect]++) effectlist[gene]=effectlist[gene] (effectlist[gene]?",":"") effect
    if (!seen_bio[gene,biotype]++) biolist[gene]=biolist[gene] (biolist[gene]?",":"") biotype
}
END {
    for (gene in count) {
        printf "%s\t%d\t%.3f\t%.3f\t%s\t%s\n",
            gene, count[gene], minfreq[gene], maxfreq[gene], biolist[gene], effectlist[gene]
    }
}' "${PREFIX}.known_sites_freq10pct_genic.tsv" | sort -k2,2nr) | awk 'BEGIN{print "GENE\tEDITING_SITES\tMIN_FREQ\tMAX_FREQ\tBIOTYPE\tEFFECTS"}{print}' > "${PREFIX}_gene_level_editing.tsv"

# 4. subset to only protein-coding genes
echo "Extracting protein-coding genes only from the gene-level summary..."
# check if the input file exists and is not empty
if [[ ! -s "${PREFIX}_gene_level_editing.tsv" ]]; then
    echo "Error: input file '${PREFIX}_gene_level_editing.tsv' not found or empty" >&2
    exit 1
fi

awk -F'\t' 'NR==1 || ($2 > 1 && $6 !~ /,/ && $5 == "protein_coding")' "${PREFIX}_gene_level_editing.tsv" > "${PREFIX}_gene_level_editing_subset_protein_coding_gt2_sites.tsv"

echo "Done."
