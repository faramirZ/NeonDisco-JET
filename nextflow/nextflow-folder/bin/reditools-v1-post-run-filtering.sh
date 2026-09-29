#!/bin/bash

# REDItools DNA-RNA Variant Filtering Pipeline (Improved Version)
# Optimized for: Shallow WGS, WES, and High-Depth WGS.
# Strategy: Treat WGS data as a tool for excluding false positives (supports novel site discovery) 
#           rather than a barrier for entry (removing known RE sites due to lack of coverage).
# Key modifications:
#   1. REDIportal annotation moved earlier to enable rescue of known sites
#   2. Known REDIportal sites bypass both dbSNP and WGS coverage filters
#   3. Novel sites still require stringent WGS validation (coverage + invariance)
#   4. Supports datasets with variable WGS coverage (0x to high depth)

set -e  
set -u  
set -o pipefail 

# ============================================================================
# CONFIGURATION
# ============================================================================

if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <RNA_BAM> <EDITING_DIR> <GENOME_FA> <SPLICE_SITES> <NOCHR_FILE> <RMSK_GTF> <SNP_GTF> <REDIPORTAL_GTF> <THREADS> <SAMPLE_ID>"
    exit 1
fi

RNA_BAM=$1              # Aligned.sortedByCoord.out.bam
EDITING_DIR=$2          # Output folder of REDItools DNA-RNA run
GENOME_FA=$3            # GRCh38.primary_assembly.genome.fa
SPLICE_SITES=$4         # gencode.v49.primary_assembly.annotation.splicesites.txt
EXClUDE_CONTIGS=$5      # excluded_contigs.txt
RMSK_GTF=$6             # rmsk.sorted.gtf.gz
SNP_GTF=$7              # snp151.sorted.gtf.gz
REDIPORTAL_GTF=$8       # REDIPortal_hg_38_atlas.gtf.gz
THREADS=$9
SAMPLE_ID=${10} 

# ============================================================================
# MAIN PIPELINE
# ============================================================================

echo "========================================================================"
echo "REDItools DNA-RNA Variant Filtering Pipeline (Improved)"
echo "========================================================================"
echo "Starting: $(date)"
echo ""

# ============================================================================
# Step 00: Validate and Copy Input Files
# ============================================================================

echo "[Step 00] Checking REDItools output directory..."
echo "Current directory: $(pwd)"
echo "Sample ID: ${SAMPLE_ID}"

# Check for outTable file
if ls "${EDITING_DIR}"/outTable_* 1> /dev/null 2>&1; then
    echo "Found outTable file(s):"
    ls "${EDITING_DIR}"/outTable_*
    NUM_FILES=$(ls "${EDITING_DIR}"/outTable_* | wc -l)
    if [ "${NUM_FILES}" -gt 1 ]; then
        echo "ERROR: Multiple outTable files found (${NUM_FILES}). This should not happen. Exiting."
        exit 1
    fi
    cp "${EDITING_DIR}"/outTable_* "outTable_${SAMPLE_ID}"
    echo "✓ Copied outTable file to outTable_${SAMPLE_ID}"
else
    echo "ERROR: No outTable file found in ${EDITING_DIR}. Exiting."
    exit 1
fi

# ============================================================================
# Step 01: Initial Filtering - RNA Variants with Invariant WGS data or Zero WGS Coverage
# ============================================================================
# Keep sites where:
#   - RNA shows substitution ($8 != "-")
#   - WGS is either invariant ($13 == "-") OR has no coverage ($10 == 0)
# This allows:
#   (a) Sites with WGS support showing no variation (any coverage depth)
#   (b) Sites with no WGS coverage (common in shallow WGS/WES) - validated later
# ============================================================================

echo "[Step 01] Filtering for RNA variants with invariant base in the WGS data or zero WGS read coverage..."
awk -v FS="\t" '{if ($8!="-" && ($13=="-" || $10==0)) print}' "outTable_${SAMPLE_ID}" > "outTable_${SAMPLE_ID}.filt.out"

NUM_INITIAL=$(wc -l < "outTable_${SAMPLE_ID}.filt.out")
echo "✓ Sites passing initial filter: ${NUM_INITIAL}"

# ============================================================================
# Step 02: Annotate with RepeatMasker
# ============================================================================

echo "[Step 02] Annotating with RepeatMasker..."
AnnotateTable.py -a "${RMSK_GTF}" -n rmsk -i "outTable_${SAMPLE_ID}.filt.out" -o "outTable_${SAMPLE_ID}.filt.out.rmsk" -u
if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk" ]; then
    echo "ERROR: RepeatMasker annotation failed. Exiting."
    exit 1
fi
echo "✓ RepeatMasker annotation complete"

# ============================================================================
# Step 03: Annotate with dbSNP
# ============================================================================

echo "[Step 03] Annotating with dbSNP..."
AnnotateTable.py -a "${SNP_GTF}" -n snp151 -i "outTable_${SAMPLE_ID}.filt.out.rmsk" -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" -u
if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" ]; then
    echo "ERROR: dbSNP annotation failed. Exiting."
    exit 1
fi
echo "✓ dbSNP annotation complete"

# ============================================================================
# Step 04: Annotate with REDIportal (CRITICAL - Done BEFORE filtering)
# ============================================================================
# This step is moved earlier compared to original pipeline to enable
# rescue of known editing sites that overlap dbSNP entries
# ============================================================================

echo "[Step 04] Annotating with REDIportal database..."

if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" ]; then
    echo "ERROR: Annotated file outTable_${SAMPLE_ID}.filt.out.rmsk.snp not found. Exiting."
    exit 1
fi

AnnotateTable.py -a "${REDIPORTAL_GTF}" -n ed -k R -c 1 -i "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed" -u
if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed" ]; then
    echo "ERROR: REDIportal annotation failed. Exiting."
    exit 1
fi

# Use this fully annotated file for all downstream steps
INPUT_FILE="outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed"
echo "✓ REDIportal annotation complete"

# ============================================================================
# Step 05: Create RNA-seq Quality Position Sets
# ============================================================================

echo "[Step 05] Selecting RNA-seq quality sets..."

# Set 1: ≥5 RNAseq reads, ≥1 mismatch, any frequency
echo "  - Set 1: ≥5 RNAseq reads, ≥1 mismatch"
selectPositions.py -i "$INPUT_FILE" -c 5 -v 1 -f 0.0 -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1"
NUM_SET1=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1")
echo "    ✓ Set 1: ${NUM_SET1} positions"

# Set 2: ≥10 RNAseq reads, ≥3 mismatches, freq ≥0.1
echo "  - Set 2: ≥10 RNAseq reads, ≥3 mismatches, freq ≥0.1"
selectPositions.py -i "$INPUT_FILE" -c 10 -v 3 -f 0.1 -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2"
NUM_SET2=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2")
echo "    ✓ Set 2: ${NUM_SET2} positions"

# ============================================================================
# Step 06-A: Select ALU Sites
# ============================================================================
# Filtering logic:
#   Known sites (REDIportal): Bypass WGS coverage and dbSNP filters
#   Novel sites: Require ≥10 WGS reads, WGS invariance, not in dbSNP
# ============================================================================

echo "[Step 06-A] Selecting ALU sites (relaxed for known, stringent for novel)..."
awk -v FS="\t" '{
    is_known = ($19=="ed");
    is_snp = ($17!="-");
    
    # WGS validation: Known sites bypass; novel sites need ≥10 reads AND invariance
    pass_wgs = (is_known || ($10>=10 && $13=="-"));
    
    # SNP filter: Known sites bypass; novel sites must not be in dbSNP
    pass_snp = (is_known || !is_snp);
    
    # Main filter: non-mitochondrial, ALU repeat, has RNA variant, passes filters
    if ($1!="chrM" && substr($16,1,3)=="Alu" && $8!="-" && pass_wgs && pass_snp) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU"

NUM_ALU=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU")
echo "✓ ALU sites selected: ${NUM_ALU}"

# ============================================================================
# Step 06-B: Select REP NON ALU Sites
# ============================================================================
# Filtering logic:
#   Known sites: Bypass WGS and dbSNP filters
#   Novel sites: Require ≥10 WGS reads, WGS invariance, low WGS freq (≤0.05), not in dbSNP
# Additional filters: Exclude simple repeats and low complexity regions
# ============================================================================

echo "[Step 06-B] Selecting REP NON ALU sites (relaxed for known, stringent for novel)..."
awk -v FS="\t" '{
    is_known = ($19=="ed");
    is_snp = ($17!="-");
    
    # WGS validation: Known sites bypass; novel sites need coverage, invariance, AND low freq
    pass_wgs = (is_known || ($10>=10 && $13=="-" && $14<=0.05));
    
    # SNP filter: Known sites bypass; novel sites must not be in dbSNP
    pass_snp = (is_known || !is_snp);
    
    # Main filter: non-mitochondrial, non-ALU repeat, excluding simple/low-complexity,
    #              has RNA variant, editing freq ≥0.1, passes validation filters
    if ($1!="chrM" && 
        substr($16,1,3)!="Alu" && 
        $15!="-" && 
        $15!="Simple_repeat" && 
        $15!="Low_complexity" && 
        $8!="-" && 
        $9>=0.1 && 
        pass_wgs && 
        pass_snp) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU"

NUM_NONALU=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU")
echo "✓ REP NON ALU sites selected: ${NUM_NONALU}"

# ============================================================================
# Step 06-C: Select NON REP Sites
# ============================================================================
# Filtering logic:
#   Known sites: Bypass WGS and dbSNP filters
#   Novel sites: Require ≥10 WGS reads, WGS invariance, low WGS freq (≤0.05), not in dbSNP
# These are sites outside annotated repetitive elements
# ============================================================================

echo "[Step 06-C] Selecting NON REP sites (relaxed for known, stringent for novel)..."
awk -v FS="\t" '{
    is_known = ($19=="ed");
    is_snp = ($17!="-");
    
    # WGS validation: Known sites bypass; novel sites need coverage, invariance, AND low freq
    pass_wgs = (is_known || ($10>=10 && $13=="-" && $14<=0.05));
    
    # SNP filter: Known sites bypass; novel sites must not be in dbSNP
    pass_snp = (is_known || !is_snp);
    
    # Main filter: non-mitochondrial, not in repeat regions, has RNA variant,
    #              editing freq ≥0.1, passes validation filters
    if ($1!="chrM" && 
        substr($16,1,3)!="Alu" && 
        $15=="-" && 
        $8!="-" && 
        $9>=0.1 && 
        pass_wgs && 
        pass_snp) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP"

NUM_NONREP=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP")
echo "✓ NON REP sites selected: ${NUM_NONREP}"

# ============================================================================
# Step 07: Extract and Consolidate Known Editing Events
# ============================================================================

echo "[Step 07] Extracting and consolidating known editing events..."

# Verify all three files exist
if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU" ] || \
   [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU" ] || \
   [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP" ]; then
    echo "ERROR: One or more filtered files are missing. Exiting."
    exit 1
fi

# Set variables
ALU_SITES="outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU"
NONALU_SITES="outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU"
NONREP_SITES="outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP"

# make dir
mkdir -p upstream-site-level-split-tables

# Combine all categories
cat "${ALU_SITES}" "${NONALU_SITES}" "${NONREP_SITES}" > "${SAMPLE_ID}-ALU-NONALU-NONREP"

# Extract known editing sites
awk -v FS="\t" '{if ($19=="ed") print}' "${SAMPLE_ID}-ALU-NONALU-NONREP" > "upstream-site-level-split-tables/${SAMPLE_ID}-knownEditing.txt"

NUM_KNOWN=$(wc -l < "upstream-site-level-split-tables/${SAMPLE_ID}-knownEditing.txt")
echo "✓ Known editing sites recovered: ${NUM_KNOWN}"

# ============================================================================
# Step 08: Process Novel Editing Candidates
# ============================================================================
# Novel sites are processed separately for:
#   A) ALU sites (less stringent re-analysis)
#   B) REP NON ALU + NON REP sites (stringent re-analysis)
# ============================================================================

echo "[Step 08] Processing novel editing candidates..."

# 08-A: Extract novel REP NON ALU and NON REP sites
echo "  [08-A] Extracting novel REP NON ALU and NON REP sites..."
cat "${NONALU_SITES}" "${NONREP_SITES}" | awk -v FS="\t" '{if ($19!="ed") print}' - > "upstream-site-level-split-tables/${SAMPLE_ID}-NONALU-NONREP--novelEditing.txt"

NUM_NOVEL_NONALU_NONREP=$(wc -l < "upstream-site-level-split-tables/${SAMPLE_ID}-NONALU-NONREP--novelEditing.txt")
echo "    ✓ Novel NONALU+NONREP sites: ${NUM_NOVEL_NONALU_NONREP}"

# Skip if no novel sites
if [ "${NUM_NOVEL_NONALU_NONREP}" -eq 0 ]; then
    echo "    ⚠ No novel NONALU+NONREP sites found. Skipping Steps 09-B through 12."
    SKIP_NONALU_NONREP=true
else
    SKIP_NONALU_NONREP=false
    # Convert to GFF
    echo "  [08-A] Converting to GFF format..."
    TableToGFF.py -i "upstream-site-level-split-tables/${SAMPLE_ID}-NONALU-NONREP--novelEditing.txt" -s -t -o "${SAMPLE_ID}-NONALU-NONREP--novelEdits.gff"
    if [ ! -f "${SAMPLE_ID}-NONALU-NONREP--novelEdits.gff" ]; then
        echo "ERROR: Failed to create GFF file for NONALU & NONREP sites. Exiting."
        exit 1
    fi
    echo "    ✓ GFF conversion complete"
fi

# 08-B: Extract novel ALU sites
echo "  [08-B] Extracting novel ALU sites..."
awk -v FS="\t" '{if ($19!="ed") print}' "${ALU_SITES}" > "upstream-site-level-split-tables/${SAMPLE_ID}-ALU--novelEditing.txt"

NUM_NOVEL_ALU=$(wc -l < "upstream-site-level-split-tables/${SAMPLE_ID}-ALU--novelEditing.txt")
echo "    ✓ Novel ALU sites: ${NUM_NOVEL_ALU}"

# Skip if no novel ALU sites
if [ "${NUM_NOVEL_ALU}" -eq 0 ]; then
    echo "    ⚠ No novel ALU sites found. Skipping Step 09-A."
    SKIP_ALU=true
else
    SKIP_ALU=false
    # Convert to GFF
    echo "  [08-B] Converting to GFF format..."
    TableToGFF.py -i "upstream-site-level-split-tables/${SAMPLE_ID}-ALU--novelEditing.txt" -s -t -o "${SAMPLE_ID}-ALU--novelEdits.gff"
    if [ ! -f "${SAMPLE_ID}-ALU--novelEdits.gff" ]; then
        echo "ERROR: Failed to create GFF file for ALU. Exiting."
        exit 1
    fi
    echo "    ✓ GFF conversion complete"
fi

# ============================================================================
# Step 09-A: Re-analyze Novel ALU Sites with Stringent Criteria
# ============================================================================
# clean up just in case 
rm -rf -- "${SAMPLE_ID}-firstalu"

if [ "${SKIP_ALU}" = false ]; then
    echo "[Step 09-A] Re-analyzing novel ALU sites with stringent criteria..."
    
    if [ ! -f "${SAMPLE_ID}-ALU--novelEdits.sorted.gff.gz" ]; then
        echo "ERROR: GFF file ${SAMPLE_ID}-ALU--novelEdits.sorted.gff.gz not found. Exiting."
        exit 1
    fi
    
    echo "  ✓ GFF file found. Running REDItoolDnaRna.py..."
    
    REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} \
        -i "${RNA_BAM}" \
        -f "${GENOME_FA}" \
        -c 5,5 -q 30,30 -m 255,255 -O 5,5 \
        -p -u -a 11-6 -l -v 1 -n 0.0 -e \
        -T "${SAMPLE_ID}-ALU--novelEdits.sorted.gff.gz" \
        -w "${SPLICE_SITES}" \
        -k "${EXClUDE_CONTIGS}" \
        -R -o "${SAMPLE_ID}-firstalu"
    
    echo "  ✓ Novel ALU re-analysis complete"
else
    echo "[Step 09-A] SKIPPED - No novel ALU sites"
fi

# ============================================================================
# Step 09-B: Re-analyze Novel NONALU+NONREP Sites and Extract Reads
# ============================================================================
# clean up just in case 
rm -rf -- "${SAMPLE_ID}-first"

if [ "${SKIP_NONALU_NONREP}" = false ]; then
    echo "[Step 09-B] Re-analyzing novel NONALU+NONREP sites with stringent criteria..."
    
    if [ ! -f "${SAMPLE_ID}-NONALU-NONREP--novelEdits.sorted.gff.gz" ]; then
        echo "ERROR: GFF file ${SAMPLE_ID}-NONALU-NONREP--novelEdits.sorted.gff.gz not found. Exiting."
        exit 1
    fi
    
    echo "  ✓ GFF file found. Running REDItoolDnaRna.py..."
    
    REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} \
        -i "${RNA_BAM}" \
        -f "${GENOME_FA}" \
        -c 10,10 -q 30,30 -m 255,255 -O 5,5 \
        -p -u -a 11-6 -l -v 3 -n 0.1 -e \
        -T "${SAMPLE_ID}-NONALU-NONREP--novelEdits.sorted.gff.gz" \
        -w "${SPLICE_SITES}" \
        -k "${EXClUDE_CONTIGS}" \
        --reads -R --addP -o "${SAMPLE_ID}-first"
    
    echo "  ✓ Novel NONALU+NONREP re-analysis complete"
else
    echo "[Step 09-B] SKIPPED - No novel NONALU & NONREP sites"
fi

# ============================================================================
# Step 10: Identify Multi-mapping Reads Using BLAT
# ============================================================================

# if [ "${SKIP_NONALU_NONREP}" = false ]; then
#     echo "[Step 10] Identifying multi-mapping reads with BLAT..."
    
#     FIRST_DNARNA=$(find "${SAMPLE_ID}-first" -type d -name "DnaRna_*" | head -n1)
    
#     if [ -z "$FIRST_DNARNA" ]; then
#         echo "ERROR: Could not find DnaRna_* directory in ${SAMPLE_ID}-first. Exiting."
#         exit 1
#     fi
    
#     echo "  ✓ Using directory: ${FIRST_DNARNA}"
    
#     pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 \
#         "${GENOME_FA}" \
#         "${FIRST_DNARNA}/outReads_$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')" \
#         "${FIRST_DNARNA}/reads.psl"
    
#     if [ ! -f "${FIRST_DNARNA}/reads.psl" ]; then
#         echo "ERROR: BLAT output file reads.psl not found in ${FIRST_DNARNA}. Exiting."
#         exit 1
#     fi
    
#     readPsl.py "${FIRST_DNARNA}/reads.psl" "${FIRST_DNARNA}/badreads.txt"
    
#     NUM_BADREADS=$(wc -l < "${FIRST_DNARNA}/badreads.txt" 2>/dev/null || echo "0")
#     echo "  ✓ Multi-mapping reads identified: ${NUM_BADREADS}"
# else
#     echo "[Step 10] SKIPPED - No novel NONALU & NONREP sites"
# fi

if [ "${SKIP_NONALU_NONREP}" = false ]; then
    echo "[Step 10] Identifying multi-mapping reads with BLAT..."
    
    FIRST_DNARNA=$(find "${SAMPLE_ID}-first" -type d -name "DnaRna_*" | head -n1)
    
    if [ -z "$FIRST_DNARNA" ]; then
        echo "ERROR: Could not find DnaRna_* directory in ${SAMPLE_ID}-first. Exiting."
        exit 1
    fi
    
    echo "  ✓ Using directory: ${FIRST_DNARNA}"
    
    QUERY_READS="${FIRST_DNARNA}/outReads_$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')"
    
    # Check if the outReads file exists and has content
    if [ -s "${QUERY_READS}" ]; then
        SKIP_BLAT=false
        pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 \
            "${GENOME_FA}" \
            "${QUERY_READS}" \
            "${FIRST_DNARNA}/reads.psl"
        
        if [ ! -f "${FIRST_DNARNA}/reads.psl" ]; then
            echo "ERROR: BLAT output file reads.psl not found in ${FIRST_DNARNA}. Exiting."
            exit 1
        fi
        
        readPsl.py "${FIRST_DNARNA}/reads.psl" "${FIRST_DNARNA}/badreads.txt"
        
        NUM_BADREADS=$(wc -l < "${FIRST_DNARNA}/badreads.txt" 2>/dev/null || echo "0")
        echo "  ✓ Multi-mapping reads identified: ${NUM_BADREADS}"
    else
        echo "  ⚠ No query reads found in ${QUERY_READS} (0 sites passed stringent criteria). Skipping BLAT filtering."
        SKIP_BLAT=true
        # Create an empty badreads file so downstream consolidation steps don't break
        touch "${FIRST_DNARNA}/badreads.txt"
    fi
else
    echo "[Step 10] SKIPPED - No novel NONALU & NONREP sites"
    SKIP_BLAT=true
fi

# ============================================================================
# Step 11: Extract and Deduplicate Reads
# ============================================================================

if [ "${SKIP_NONALU_NONREP}" = false ] && [ "${SKIP_BLAT}" = false ]; then
    echo "[Step 11] Extracting and deduplicating reads..."
    
    sort -k1,1 -k2,2n -k3,3n "${FIRST_DNARNA}/outPosReads_$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')" | \
        mergeBed > bed
    
    echo "  - Extracting reads from regions..."
    samtools view -@ ${THREADS} -L bed -h -b "${RNA_BAM}" > "${SAMPLE_ID}_bed.bam"
    
    echo "  - Sorting by name for fixmate..."
    samtools sort -@ ${THREADS} -n "${SAMPLE_ID}_bed.bam" -o "${SAMPLE_ID}_bed_ns.bam"
    
    echo "  - Fixing mate information..."
    samtools fixmate -@ ${THREADS} -m "${SAMPLE_ID}_bed_ns.bam" "${SAMPLE_ID}_bed_ns_fx.bam"
    
    echo "  - Sorting by coordinate..."
    samtools sort -@ ${THREADS} "${SAMPLE_ID}_bed_ns_fx.bam" -o "${SAMPLE_ID}_bed_ns_fx_st.bam"
    
    echo "  - Marking and removing duplicates..."
    samtools markdup -r -@ ${THREADS} "${SAMPLE_ID}_bed_ns_fx_st.bam" "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam"
    
    echo "  - Indexing deduplicated BAM..."
    samtools index "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam"
    
    echo "  ✓ Deduplication complete"
else
    echo "[Step 11] SKIPPED - No novel NONALU & NONREP sites"
fi

# ============================================================================
# Step 12: Final Re-analysis with Deduplicated Reads and Multi-mapping Filter
# ============================================================================
# clean up just in case 
rm -rf -- "${SAMPLE_ID}-second"

if [ "${SKIP_NONALU_NONREP}" = false ] && [ "${SKIP_BLAT}" = false ]; then
    echo "[Step 12] Final re-analysis with deduplicated reads and multi-mapping filter..."
    
    if [ ! -f "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam" ]; then
        echo "ERROR: Deduplicated BAM file ${SAMPLE_ID}_bed_ns_fx_st_dedup.bam not found. Exiting."
        exit 1
    fi
    
    echo "  ✓ Deduplicated BAM found. Running final REDItoolDnaRna.py..."
    
    REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} \
        -i "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam" \
        -f "${GENOME_FA}" \
        -c 10,10 -q 30,30 -m 255,255 -O 5,5 \
        -p -u -a 11-6 -l -v 3 -n 0.1 -e \
        -T "${SAMPLE_ID}-NONALU-NONREP--novelEdits.sorted.gff.gz" \
        -w "${SPLICE_SITES}" \
        -k "${EXClUDE_CONTIGS}" \
        -b "${FIRST_DNARNA}/badreads.txt" \
        -R --rmIndels \
        -o "${SAMPLE_ID}-second"
    
    echo "  ✓ Final re-analysis complete"
else
    echo "[Step 12] SKIPPED - No novel NONALU & NONREP sites"
fi

# ============================================================================
# Step 13: Consolidate Final Results with Status Labels
# ============================================================================

echo "[Step 13] Consolidating final results with editing status labels..."

# Extract novel ALU results if they exist
if [ "${SKIP_ALU}" = false ]; then
    echo "  - Extracting novel ALU sites from REDItools output..."
    FIRSTALU_DNARNA=$(find "${SAMPLE_ID}-firstalu" -type d -name "DnaRna_*" | head -n1)

    if [ -z "$FIRSTALU_DNARNA" ] || ! ls "${FIRSTALU_DNARNA}"/outTable_* 1> /dev/null 2>&1; then
        echo "ERROR: ALU re-analysis ran but produced no outTable for sample ${SAMPLE_ID}." >&2
        exit 1
    fi

    awk 'BEGIN {OFS="\t"} !/^Region/ && !/^#/' "${FIRSTALU_DNARNA}"/outTable_* > "${SAMPLE_ID}-novelEditing-ALU.raw"

    NUM_NOVEL_ALU_FINAL=$(wc -l < "${SAMPLE_ID}-novelEditing-ALU.raw")

else
    NUM_NOVEL_ALU_FINAL=0
fi

# Extract novel NONALU+NONREP results if they exist
if [ "${SKIP_NONALU_NONREP}" = false ] && [ "${SKIP_BLAT}" = false ]; then
    echo "  - Extracting novel NONALU+NONREP sites from REDItools output..."
    SECOND_DNARNA=$(find "${SAMPLE_ID}-second" -type d -name "DnaRna_*" | head -n1)
    
    if [ -z "$SECOND_DNARNA" ] || ! ls "${SECOND_DNARNA}"/outTable_* 1> /dev/null 2>&1; then
        echo "ERROR: NONALU+NONREP re-analysis ran but produced no outTable for sample ${SAMPLE_ID}." >&2
        exit 1
    fi

    # Copy raw output, skip any header lines
    awk 'BEGIN {OFS="\t"} !/^Region/ && !/^#/' "${SECOND_DNARNA}"/outTable_* > "${SAMPLE_ID}-novelEditing-NONALU-NONREP.raw"
    NUM_NOVEL_NONALU_NONREP_FINAL=$(wc -l < "${SAMPLE_ID}-novelEditing-NONALU-NONREP.raw")
    echo "    ✓ Novel NONALU+NONREP sites extracted: ${NUM_NOVEL_NONALU_NONREP_FINAL}"
else
    echo "  ⚠ Skipping novel NONALU+NONREP extraction (no sites survived stringent filters)."
    NUM_NOVEL_NONALU_NONREP_FINAL=0
fi

# Combine all novel sites for annotation
echo "  - Combining all novel editing sites..."
if [ "${SKIP_ALU}" = false ] && [ "${SKIP_NONALU_NONREP}" = false ] && [ "${SKIP_BLAT}" = false ]; then
    cat "${SAMPLE_ID}-novelEditing-ALU.raw" "${SAMPLE_ID}-novelEditing-NONALU-NONREP.raw" > "${SAMPLE_ID}-novelEditing.raw"
elif [ "${SKIP_ALU}" = false ]; then
    cp "${SAMPLE_ID}-novelEditing-ALU.raw" "${SAMPLE_ID}-novelEditing.raw"
elif [ "${SKIP_NONALU_NONREP}" = false ] && [ "${SKIP_BLAT}" = false ]; then
    cp "${SAMPLE_ID}-novelEditing-NONALU-NONREP.raw" "${SAMPLE_ID}-novelEditing.raw"
else
    touch "${SAMPLE_ID}-novelEditing.raw"  # Safely creates an empty file if nothing survived
fi

NUM_NOVEL_TOTAL=$(wc -l < "${SAMPLE_ID}-novelEditing.raw" 2>/dev/null || echo "0")
echo "    ✓ Total novel editing sites (before annotation): ${NUM_NOVEL_TOTAL}"

# Annotate novel sites to match known sites column structure
if [ ${NUM_NOVEL_TOTAL} -gt 0 ]; then
    echo "  - Annotating novel sites with RepeatMasker..."
    AnnotateTable.py -a "${RMSK_GTF}" -n rmsk -i "${SAMPLE_ID}-novelEditing.raw" -o "${SAMPLE_ID}-novelEditing.rmsk" -u
    
    echo "  - Annotating novel sites with dbSNP..."
    AnnotateTable.py -a "${SNP_GTF}" -n snp151 -i "${SAMPLE_ID}-novelEditing.rmsk" -o "${SAMPLE_ID}-novelEditing.rmsk.snp" -u
    
    echo "  - Adding REDIportal placeholder column (all novel sites = '-')..."
    awk 'BEGIN {OFS="\t"} {print $0, "-"}' "${SAMPLE_ID}-novelEditing.rmsk.snp" > "${SAMPLE_ID}-novelEditing.annotated"
    
    echo "    ✓ Novel sites annotation complete"
else
    echo "    ⚠ No novel sites to annotate"
    touch "${SAMPLE_ID}-novelEditing.annotated"
fi

# Add editing status labels to novel sites based on their category
echo "  - Adding editing status labels to novel sites..."
mkdir -p final-labeled-calls

if [ ${NUM_NOVEL_TOTAL} -gt 0 ]; then
    awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        # 1. Default fallback
        status = "UNKNOWN" 
        # 2. Strict Alu (Must be SINE AND start with Alu)
        if ($15=="SINE" && substr($16,1,3)=="Alu") {
            status = "NOVEL_ALU"
        }  
        # 3. Non-Rep (Must be exactly - and -)
        else if ($15=="-" && $16=="-") {
            status = "NOVEL_NONREP"
        } 
        # 4. Everything else that is a repeat but not a SINE/Alu
        else {
            status = "NOVEL_NONALU"
        }
        print $0, status
    }' "${SAMPLE_ID}-novelEditing.annotated" > final-labeled-calls/"${SAMPLE_ID}-novelEditing.tsv"

    echo "    ✓ Status labels added"
else
    touch final-labeled-calls/"${SAMPLE_ID}-novelEditing.tsv"
fi

# Process known editing sites - add status label based on category
echo "  - Labeling known editing sites by category..."
awk 'BEGIN {FS="\t"; OFS="\t"} 
{
    # 1. Default fallback
    status = "UNKNOWN" 

    # 2. Strict Alu (Must be SINE AND start with Alu)
    if ($15=="SINE" && substr($16,1,3)=="Alu") {
        status = "KNOWN_ALU"
    } 
    # 3. Non-Rep (Must be exactly - and -)
    else if ($15=="-" && $16=="-") {
        status = "KNOWN_NONREP"
    } 
    # 4. Everything else that is a repeat but not a SINE/Alu
    else {
        status = "KNOWN_NONALU"
    }
    print $0, status
}' "upstream-site-level-split-tables/${SAMPLE_ID}-knownEditing.txt" > final-labeled-calls/"${SAMPLE_ID}-knownEditing-labeled.tsv"

# Create header-only file with EditingStatus column
echo "  - Creating final file with header..."
echo -e "Region\tPosition\tReference\tStrand\tCoverage-q\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\tRepeatType\tRepeatName\tSNPFlag\tdbSNP_ID\tREDIPortalKnownEditingSites\tEditingStatus" > final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv"

# Combine known and novel sites
cat final-labeled-calls/"${SAMPLE_ID}-knownEditing-labeled.tsv" final-labeled-calls/"${SAMPLE_ID}-novelEditing.tsv" >> final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv"

NUM_ALL_EDITING=$((NUM_KNOWN + NUM_NOVEL_TOTAL))
echo "    ✓ Total editing sites (known + novel): ${NUM_ALL_EDITING}"

# Create a summary table
echo "  - Generating summary statistics..."
mkdir -p summary-stats
echo -e "EditingStatus\tCount" > summary-stats/"${SAMPLE_ID}-editingSummary.txt"
awk -F"\t" 'NR>1 {print $NF}' final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv" | sort | uniq -c | \
    awk '{print $2"\t"$1}' >> summary-stats/"${SAMPLE_ID}-editingSummary.txt"

echo "✓ Result consolidation complete"
echo ""
echo "Final output file: ${SAMPLE_ID}-allEditing (with header)"
echo "Summary statistics: ${SAMPLE_ID}-editingSummary.txt"

# ============================================================================
# Step 14: Final Filter, Category Discovery & Deep Statistics
# ============================================================================
# echo "[Step 14] Generating AG-Substitution-Site-only list and comprehensive statistics..."

# # 14-A: Create the High-Confidence Master List
# # Filter: Canonical A-to-G only
# awk 'BEGIN {FS="\t"; OFS="\t"} 
# NR==1 {print $0} 
# NR>1 {
#     if ($3=="A" && $8=="AG") print $0
# }' final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv" > final-labeled-calls/"${SAMPLE_ID}-AG-Subs-Only-Sites.tsv"

# # 14-B: Site-Based Category Discovery
# # This counts unique genomic locations using the custom labels
# TOTAL_SITES=$(grep -v "Region" final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv" | wc -l)
# TOTAL_AG=$(grep -v "Region" final-labeled-calls/"${SAMPLE_ID}-AG-Subs-Only-Sites.tsv" | wc -l)

# # Categorize specifically for the log output
# ALU_AG=$(grep -P "\t[^\t]*ALU[^\t]*$" final-labeled-calls/"${SAMPLE_ID}-AG-Subs-Only-Sites.tsv" | wc -l)
# NONALU_AG=$(grep -P "\t[^\t]*NONALU[^\t]*$" final-labeled-calls/"${SAMPLE_ID}-AG-Subs-Only-Sites.tsv" | wc -l)
# NONREP_AG=$(grep -P "\t[^\t]*NONREP[^\t]*$" final-labeled-calls/"${SAMPLE_ID}-AG-Subs-Only-Sites.tsv" | wc -l)
# AG_RATIO=$(awk -v t="$TOTAL_SITES" -v a="$TOTAL_AG" 'BEGIN {if (t > 0) printf "%.2f", (a/t)*100; else print "0"}')

# # 14-C: Run the REDItools Weighted Statistics Utility
# # This provides the deep read-level substitution distribution math
# echo "  - Running getStatistics.py for weighted read distribution..."
# get_statistics-v2.py --input final-labeled-calls/"${SAMPLE_ID}-allEditing.tsv" --output final-labeled-calls/"${SAMPLE_ID}-editingStats.txt"

# # 14-D: Generate Consolidated Console Report
# {
#     printf "%s\n" "===================================================="
#     printf "%s\n" "       RNA EDITING SITE DISCOVERY SUMMARY           "
#     printf "%s\n" "===================================================="
#     printf "Sample ID:             %s\n" "${SAMPLE_ID}"
#     printf "Total Candidate Sites: %s\n" "${TOTAL_SITES}"
#     printf "High-Conf AG Sites:    %s (Freq >= 0.1)\n" "${TOTAL_AG}"
#     printf "AG Specificity Ratio:  %s%%\n" "${AG_RATIO}"
#     printf "%s\n" "----------------------------------------------------"
#     printf "%s\n" "AG Site Distribution (Unique Locations):"
#     printf "  - Alu Elements:      %s\n" "${ALU_AG}"
#     printf "  - Non-Alu Repeats:   %s\n" "${NONALU_AG}"
#     printf "  - Non-Rep Regions:   %s\n" "${NONREP_AG}"
#     printf "%s\n" "----------------------------------------------------"
#     printf "Weighted Read Stats:   %s-editingStats.txt\n" "${SAMPLE_ID}"
#     printf "%s\n" "===================================================="
# } > "summary-stats/${SAMPLE_ID}-FinalDiscoverySummary.txt"

# safety scrub:
# - s/\x1b//g         : Removes the 'ESC' character itself
# - s/\x1b([A-Z]//g   : Removes set transitions like (B or (O
# - s/\[[0-9;]*m//g   : Removes standard color codes
# - s/[[:cntrl:]]//g  : Removes any other non-printable control characters
# sed -i 's/\x1b([A-Z]//g; s/\x1b\[[0-9;]*[a-zA-Z]//g; s/[[:cntrl:]]//g' summary-stats/"${SAMPLE_ID}-FinalDiscoverySummary.txt"

# cat summary-stats/"${SAMPLE_ID}-FinalDiscoverySummary.txt"
# echo "✓ Summary generated: ${SAMPLE_ID}-FinalDiscoverySummary.txt"
# echo "✓ Master AG list for SnpEff: ${SAMPLE_ID}-AG-Subs-Only-Sites.tsv"

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

echo ""
echo "========================================================================"
echo "Pipeline Completed Successfully!"
echo "========================================================================"
echo "Finished: $(date)"
echo ""
echo "SUMMARY OF IMPROVEMENTS FROM ORIGINAL PIPELINE:"
echo "------------------------------------------------"
echo "  1. Step 01: Relaxed WGS coverage requirement (supports shallow WGS/WES)"
echo "     - Allows sites with 0 WGS coverage to pass initial filter"
echo "     - Known sites validated by REDIportal regardless of WGS depth"
echo ""
echo "  2. Step 04: REDIportal annotation moved earlier in workflow"
echo "     - Known sites annotated BEFORE filtering steps"
echo "     - Prevents loss of bona fide editing sites overlapping dbSNP"
echo ""
echo "  3. Steps 06A-C: Differential filtering for known vs. novel sites"
echo "     - KNOWN SITES: Bypass WGS coverage AND dbSNP filters"
echo "     - NOVEL SITES: Require stringent WGS validation"
echo ""
echo "RESULTS:"
echo "--------"
echo "  Known editing sites: ${NUM_KNOWN}"
echo "  Novel ALU sites: ${NUM_NOVEL_ALU_FINAL}"
echo "  Novel NONALU+NONREP sites: ${NUM_NOVEL_NONALU_NONREP_FINAL}"
echo "  Total novel sites: ${NUM_NOVEL_TOTAL}"
echo "  TOTAL EDITING SITES: ${NUM_ALL_EDITING}"
echo ""
echo "  Category breakdown:"
echo "    Total ALU sites (known + novel): ${NUM_ALU}"
echo "    Total NONALU sites: ${NUM_NONALU}"
echo "    Total NONREP sites: ${NUM_NONREP}"
echo ""
echo "OUTPUT FILES:"
echo "-------------"
echo "  ★ FINAL COMPLETE SET (with header): ${SAMPLE_ID}-allEditing.tsv"
echo ""
echo "  Intermediate files:"
echo "    Known editing sites: ${SAMPLE_ID}-knownEditing-labeled.tsv"
echo "    Novel editing sites: ${SAMPLE_ID}-novelEditing.tsv"
echo "    Novel ALU sites: ${SAMPLE_ID}-novelEditing-ALU.raw"
echo "    Novel NONALU+NONREP sites: ${SAMPLE_ID}-novelEditing-NONALU-NONREP.raw"
echo ""
echo "========================================================================"



