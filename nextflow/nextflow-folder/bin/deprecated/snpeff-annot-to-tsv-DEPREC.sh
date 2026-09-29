#!/usr/bin/env bash
# snpeff-annot-to-tsv.sh
# Extract and flatten SnpEff-annotated transcript variants to a Transcript-Level summary matrix.
# Computes a coverage-weighted mean editing frequency for global editing-loss analysis
# and preserves maximum frequency and impact metrics for downstream analytical forks.
set -euo pipefail

# Check for the 2 required positional parameters
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.ann.vcf sample_name" >&2
    exit 1
fi

INPUT="$1"
SAMPLE_NAME="$2"

# Validation
if [[ ! -f "$INPUT" ]]; then
    echo "Error: input VCF file '$INPUT' not found" >&2
    exit 1
fi

# Define the explicit static column mapping layout for SnpSift extractFields
# Field positions inside awk:
# $1=CHROM, $2=POS, $3=REF, $4=ALT, $5=ANN, $6=STRAND, $7=FREQ, $8=COV, $9=REPTYPE, $10=REPNAME, $11=EDITSTATUS
FIELDS="CHROM POS REF ALT ANN STRAND FREQ COV REPTYPE REPNAME EDITSTATUS"

# -------------------------------------------------------------------------
# Step 1: All known editing sites (Master Extraction Pass)
# Write to a LOCAL TEMPORARY file first
# -------------------------------------------------------------------------

echo "Extracting all known sites..."
SnpSift filter "REDIPORTAL == 'ed'" "$INPUT" | \
SnpSift extractFields - $FIELDS \
> "${SAMPLE_NAME}.known_sites_all_freq.tsv.part"

# Validate the working file
if [[ ! -s "${SAMPLE_NAME}.known_sites_all_freq.tsv.part" ]]; then
    echo "Error: The content of ${SAMPLE_NAME} appears to be empty" >&2
    rm -f "${SAMPLE_NAME}.known_sites_all_freq.tsv.part"
    exit 1
fi

# Only now does the "real" output name come into existence
mv "${SAMPLE_NAME}.known_sites_all_freq.tsv.part" \
   "${SAMPLE_NAME}.known_sites_all_freq.tsv"

# -------------------------------------------------------------------------
# Step 2: Collapse to Transcript-Level Matrices via In-Memory Stream Processing
# -------------------------------------------------------------------------
echo "Collapsing transcript-level annotations and calculating weights..."

# Map positional indices: 1=TX, 2=GENE, 3=SITES, 4=WEIGHTED_MEAN, 5=MAX_FREQ, 6=BIOTYPE, 7=IMPACT, 8=EFFECTS
awk -F'\t' '
NR>1 {
    site_id = $1 ":" $2
    freq    = $7 + 0.0
    cov     = $8 + 0.0
    split($5, transcripts, ",")
    for (t in transcripts) {
        split(transcripts[t], subfields, "|")
        effect  = subfields[2]
        impact  = subfields[3]   # Impact (HIGH, MODERATE, LOW, MODIFIER)
        gene    = subfields[4]
        tx_id   = subfields[7]
        biotype = subfields[8]
        if (tx_id == "" || tx_id == "." || effect == "intergenic_region" || biotype == "intergenic_region") {
            continue
        }
        if (!seen_site[tx_id, site_id]++) {
            count[tx_id]++
        }
        gene_name[tx_id] = gene
        biolist[tx_id]   = biotype
        
        sum_freq_cov[tx_id] += (freq * cov)
        sum_cov[tx_id]      += cov
        
        if (freq > maxfreq[tx_id]) {
            maxfreq[tx_id] = freq
        }
        if (effect != "" && !seen_effect[tx_id, effect]++) {
            effectlist[tx_id] = (effectlist[tx_id] ? effectlist[tx_id] "," : "") effect
        }
        if (impact != "" && !seen_impact[tx_id, impact]++) {
            impactlist[tx_id] = (impactlist[tx_id] ? impactlist[tx_id] "," : "") impact
        }
    }
}
END {
    for (tx in count) {
        weighted_mean = (sum_cov[tx] > 0) ? (sum_freq_cov[tx] / sum_cov[tx]) : 0.0
        printf "%s\t%s\t%d\t%.4f\t%.4f\t%s\t%s\t%s\n",
            tx, gene_name[tx], count[tx], weighted_mean, maxfreq[tx], biolist[tx], impactlist[tx], effectlist[tx]
    }
}' "${SAMPLE_NAME}.known_sites_all_freq.tsv" | sort -k3,3nr -k4,4nr > "${SAMPLE_NAME}_transcript_level_editing.tmp"

# -------------------------------------------------------------------------
# Step 3: Downstream Analytical Splits and Specialized Sorting Forks
# Safe Atomic Generation (Rename/Move at the very last second)
# -------------------------------------------------------------------------

awk 'BEGIN{print "TRANSCRIPT\tGENE\tEDITING_SITES\tWEIGHTED_MEAN_FREQ\tMAX_FREQ\tBIOTYPE\tIMPACT\tEFFECTS"}{print}' \
    "${SAMPLE_NAME}_transcript_level_editing.tmp" \
    > "${SAMPLE_NAME}_transcript_level_editing.tsv"

# User-configurable filtering threshold for downstream neoantigen execution
NEO_MIN_FREQ=0.2000

# --- FORK 1: Hyper-Editing (Clustered sites, sorted by global burden) ---
# Filters for transcripts with high-density clusters of editing sites (>2) to capture structural/stability shifts.
awk -F'\t' '
    NR == 1 { next }
    $6 == "protein_coding" && $3 > 2 { print }
' "${SAMPLE_NAME}_transcript_level_editing.tsv" | \
sort -t$'\t' -k4,4nr > "${SAMPLE_NAME}_hyper_edited_protein-coding.tmp"

awk 'BEGIN{print "TRANSCRIPT\tGENE\tEDITING_SITES\tWEIGHTED_MEAN_FREQ\tMAX_FREQ\tBIOTYPE\tIMPACT\tEFFECTS"}{print}' \
    "${SAMPLE_NAME}_hyper_edited_protein-coding.tmp" \
    > "${SAMPLE_NAME}_hyper_edited_protein-coding.tsv"


# --- FORK 2: Isolated High-Penetrance Sites (Clean, highly efficient singletons/doubletons) ---
# Isolates low-site-count protein-coding transcripts (<=2) where the single or double sites 
# are edited with high penetrance (>=50%), highlighting major localized functional or regulatory switch-flips.
awk -F'\t' '
    NR == 1 { next }
    $6 == "protein_coding" && $3 <= 2 && $5 >= 0.5000 { print }
' "${SAMPLE_NAME}_transcript_level_editing.tsv" | \
sort -t$'\t' -k5,5nr > "${SAMPLE_NAME}_isolated_high_penetrance_sites.tmp"

awk 'BEGIN{print "TRANSCRIPT\tGENE\tEDITING_SITES\tWEIGHTED_MEAN_FREQ\tMAX_FREQ\tBIOTYPE\tIMPACT\tEFFECTS"}{print}' \
    "${SAMPLE_NAME}_isolated_high_penetrance_sites.tmp" \
    > "${SAMPLE_NAME}_isolated_high_penetrance_sites_protein-coding.tsv"


# --- FORK 3: Neoantigen Candidates (Altered coding framework with frequency filtering) ---
# Isolates protein-coding transcripts carrying translational structural/amino-acid modifications (HIGH or MODERATE),
# filtered by a user-level frequency threshold to ensure high-confidence expression and clonal presentation.
awk -F'\t' -v min_freq="$NEO_MIN_FREQ" '
    NR == 1 { next }
    $6 == "protein_coding" && ($7 ~ /MODERATE/ || $7 ~ /HIGH/) && $5 >= min_freq { print }
' "${SAMPLE_NAME}_transcript_level_editing.tsv" | \
sort -t$'\t' -k5,5nr > "${SAMPLE_NAME}_neoantigen_candidates.tmp"

awk 'BEGIN{print "TRANSCRIPT\tGENE\tEDITING_SITES\tWEIGHTED_MEAN_FREQ\tMAX_FREQ\tBIOTYPE\tIMPACT\tEFFECTS"}{print}' \
    "${SAMPLE_NAME}_neoantigen_candidates.tmp" \
    > "${SAMPLE_NAME}_neoantigen_candidates_protein-coding.tsv"


# Clean up temporary scraps safely
rm -f *.tmp

echo "All files successfully generated. Pipeline complete."

