# NeonDisco-TEAL — JET Identification Pipeline V.1.1.0

```
  ███╗   ██╗███████╗ ██████╗ ███╗   ██╗██████╗ ██╗███████╗ ██████╗ ██████╗
  ████╗  ██║██╔════╝██╔═══██╗████╗  ██║██╔══██╗██║██╔════╝██╔════╝██╔═══██╗
  ██╔██╗ ██║█████╗  ██║   ██║██╔██╗ ██║██║  ██║██║███████╗██║     ██║   ██║
  ██║╚██╗██║██╔══╝  ██║   ██║██║╚██╗██║██║  ██║██║╚════██║██║     ██║   ██║
  ██║ ╚████║███████╗╚██████╔╝██║ ╚████║██████╔╝██║███████║╚██████╗╚██████╔╝
  ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚═════╝ ╚═╝╚══════╝ ╚═════╝ ╚═════╝

  ████████╗███████╗ █████╗ ██╗
  ╚══██╔══╝██╔════╝██╔══██╗██║
     ██║   █████╗  ███████║██║
     ██║   ██╔══╝  ██╔══██║██║
     ██║   ███████╗██║  ██║███████╗
     ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
```

> **JET identification and MHC-I binding prediction, adapted for human GRCh38**
> Based on Burbage & Rocañín-Arjó et al., *Science Immunology* 2023

---

## Table of Contents

1. [Background](#1-background)
2. [Pipeline Overview](#2-pipeline-overview)
3. [Repository Structure](#3-repository-structure)
4. [Software Requirements](#4-software-requirements)
5. [Reference Files Required](#5-reference-files-required)
6. [Input Data Preparation](#6-input-data-preparation)
7. [Pre-Run Checklist](#7-pre-run-checklist)
8. [Running the Pipeline](#8-running-the-pipeline)
9. [Expected Outputs](#9-expected-outputs)
10. [Interpreting the Results](#10-interpreting-the-results)
11. [Integrity Check](#11-integrity-check)
12. [Moving to Another Machine](#12-moving-to-another-machine)
13. [Troubleshooting](#13-troubleshooting)
14. [Citation](#14-citation)
15. [Authors and Contact](#15-authors-and-contact)

---

## 1. Background

Transposable elements (TEs) make up nearly half of the human genome and are normally silenced by epigenetic mechanisms such as H3K9me3 deposition via SETDB1. In tumour cells this silencing is frequently disrupted, allowing TEs to be aberrantly transcribed and spliced into gene transcripts.

**JETs (Junctions between Exons and Transposable Elements)** are non-canonical splice junctions that arise when an exonic splice site is paired with a cryptic splice site inside a TE. The resulting fusion transcript can encode a chimeric peptide — part gene product, part transposable element — that is absent from normal tissue. If such a peptide is presented on MHC class I at the tumour cell surface, it is a candidate tumour-specific neoantigen.

This pipeline identifies JETs from RNA-seq data, translates the junctions into fusion peptides (8–11 amino acids), and predicts which peptides are likely to be presented by the patient's HLA class I alleles. The JET-calling logic was developed by Alexandre Houy and Christel Goudot at Institut Curie and described in:

> Burbage M, Rocañín-Arjó A, et al. *Epigenetically-controlled tumor antigens derived from splice junctions between exons and transposable elements.* **Science Immunology**, 2023.

This repository adapts the original pipeline for **human GRCh38**, replaces the original netMHCpan step with **pVACbind** (pVACtools) running several predictors through an Apptainer image, and wraps everything in one master script (`run_TEAL_pipeline.sh`) that processes many samples per run.

---

## 2. Pipeline Overview

```
RNA-seq FASTQ (.fastq or .fastq.gz)
      │
      ▼
┌──────────────────────────────────────────────┐
│  STEP 1 — STAR 2-pass alignment              │
│  Builds the genome index if it is missing;   │
│  reads .gz directly (no decompression);      │
│  produces sorted BAM + chimeric junctions    │
└─────────────────────┬────────────────────────┘
                      │  Chimeric.out.junction, SJ.out.tab, Log.final.out
                      ▼
┌──────────────────────────────────────────────┐
│  STEP 2 — JET identification (R)             │
│  Annotates junctions against RepeatMasker,   │
│  exons and promoters; keeps repeat–exon      │
│  junctions; translates them in 3 frames      │
│  for peptide sizes 8, 9, 10 and 11           │
└─────────────────────┬────────────────────────┘
                      │  Fusions ... size{8-11}.fasta / .ids.txt
                      ▼
┌──────────────────────────────────────────────┐
│  STEP 3 — pVACbind (Apptainer image)         │
│  NetMHCpan, NetMHCpanEL, MHCflurryEL,        │
│  BigMHC_IM, PRIME, DeepImmuno (+ NetMHCstab) │
│  against the sample's HLA alleles;           │
│  one run per peptide size                    │
└─────────────────────┬────────────────────────┘
                      │  pvacbind_size{8-11}/MHC_Class_I/*.tsv
                      ▼
┌──────────────────────────────────────────────┐
│  STEP 4 — Compile results (per sample)       │
│  Combines filtered.tsv of sizes 8–11, adds   │
│  Sample, Size and Junction columns           │
└─────────────────────┬────────────────────────┘
                      ▼
        <sample>_final_combined_data_TEAL.tsv
```

For cohorts where HLA types are not yet known, run **Step 1 only**, type HLA from the aligned BAMs, then resume at Step 2 (see [8.4](#84-workflow-for-a-new-cohort-no-hla-typing-yet)).

---

## 3. Repository Structure

```
neondisco-jet/
├── JET_identification_pipeline/
│   ├── run_TEAL_pipeline.sh            # Master orchestrator (entry point)
│   ├── step1_alignment.sh              # STAR alignment (+ index build if missing)
│   ├── step2_jet_identification.sh     # Calls the R script for sizes 8–11
│   ├── JET_analysis_filtered.R         # JET annotation and peptide extraction
│   ├── step3_pvacbind.sh               # pVACbind binding prediction (via SIF)
│   ├── step4_compiler.sh               # Per-sample compiled results table
│   ├── check_integrity.sh              # md5sum comparison of two runs
│   └── lib/
│       └── common.sh                   # Shared logging / helper functions
├── infoDir/
│   └── samples.tsv                     # Sample metadata (5 columns)
├── inputData/
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   └── Homo_sapiens.GRCh38.115.gtf
├── metadataDir/
│   ├── hg38.repeatmasker_reformat.txt  # RepeatMasker table (see 5.3)
│   └── GRCh38_115.txdb                 # Cached TxDb (auto-built on first Step 2)
└── renv/                               # R package library
```

> `step3_netmhcpan.sh` and `aggregate_binders.sh` are the superseded netMHCpan-era scripts. The master script no longer calls them.

---

## 4. Software Requirements

| Tool | Notes |
|---|---|
| STAR, samtools | Installed through pixi (`pixi.toml`). The original pipeline used STAR 2.5.3a / samtools 1.3. Samtools is only used to index the aligned BAM. |
| R 4.5.x | Installed system-wide (rig). Packages live in the project `renv` library. |
| Apptainer | Runs the pVACtools image. |
| pVACtools image | `pVACtools.sif` (developed with pVACtools 7.1.2). Must contain the local IEDB install at `/opt/iedb` — this is what keeps NetMHCpanEL from falling back to the slow IEDB web API. |

The standalone `netMHCpan` installation is **no longer required**; NetMHCpan is run from inside the image.

### R packages

```r
# In R, from ~/repos/neondisco-jet/
BiocManager::install(c(
    "GenomicAlignments", "GenomicFeatures", "GenomeInfoDb",
    "txdbmaker", "AnnotationDbi",
    "biovizBase", "ggbio",
    "BSgenome.Hsapiens.UCSC.hg38"
))
install.packages(c("ggplot2", "data.table", "scales", "getopt"))
```

> `biovizBase` needs the system library `libuv1-dev`. If installation fails, ask your sysadmin to install it.

### pixi environment

```bash
cd ~/repos/neondisco-jet
pixi install
```

### pVACtools image

Place the image somewhere readable (e.g. `~/images/pVACtools.sif`). If you build it yourself, for example with `apptainer pull pVACtools.sif docker://griffithlab/pvactools:<tag>`, make sure it passes the checks in [Section 7](#7-pre-run-checklist).

---

## 5. Reference Files Required

### 5.1 Genome FASTA

```bash
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/repos/neondisco-jet/inputData/
```

### 5.2 Gene annotation GTF

```bash
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
mv Homo_sapiens.GRCh38.115.gtf ~/repos/neondisco-jet/inputData/
```

### 5.3 RepeatMasker annotation (UCSC Table Browser)

Step 2 uses this table to decide which junction ends fall in a repeat element. It is a **genome-level** file: download it once per genome build and point every sample at it. The file name is free; existing names such as `379T.repeatmasker_reformat.txt` work for any hg38 sample.

1. Go to the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
2. Set:

| Field | Value |
|---|---|
| Clade | Mammal |
| Genome | Human |
| Assembly | GRCh38/hg38 |
| Group | Variation and Repeats |
| Track | RepeatMasker |
| Table | rmsk |
| Region | genome |
| Output format | selected fields from primary and related table |
| Output file | `hg38.repeatmasker_reformat.txt` |
| File type returned | plain text |

3. Click **Get output**, then select these columns, in this order: `genoName`, `genoStart`, `genoEnd`, `strand`, `repName`, `repClass`, `repFamily`
4. Click **Get output**, save the file (keep the header line) and move it to `~/repos/neondisco-jet/metadataDir/`.

Chromosome names (`chr1` vs `1`) are reconciled automatically in the R script.

### 5.4 STAR genome index

You do not build this by hand. `step1_alignment.sh` checks for `genomeParameters.txt` inside `--star-index`; if it is missing the index is built (about 40 minutes and ~30 GB RAM, once), otherwise the build is skipped. The index directory is created automatically.

---

## 6. Input Data Preparation

### 6.1 `samples.tsv`

One row per sample, **5 tab-separated columns**, no header:

```
name<TAB>R1_path<TAB>R2_path<TAB>HLA_alleles<TAB>repeatmasker_file
```

| Column | Description | Example |
|---|---|---|
| `name` | Unique sample identifier (becomes the output folder prefix) | `sampleA` |
| `R1_path` | Absolute path to R1 FASTQ (`.fastq` or `.fastq.gz`) | `/data/raw/sampleA_R1.fastq.gz` |
| `R2_path` | Absolute path to R2 FASTQ | `/data/raw/sampleA_R2.fastq.gz` |
| `HLA_alleles` | Comma-separated class I alleles, e.g. `HLA-A*02:01,HLA-B*07:02` | see below |
| `repeatmasker_file` | Path to the RepeatMasker table (5.3) | `/home/faramir/repos/neondisco-jet/metadataDir/hg38.repeatmasker_reformat.txt` |

```
sampleA	/data/raw/sampleA_R1.fastq.gz	/data/raw/sampleA_R2.fastq.gz	HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*40:01,HLA-C*03:04,HLA-C*07:02	/home/faramir/repos/neondisco-jet/metadataDir/hg38.repeatmasker_reformat.txt
sampleB	/data/raw/sampleB_R1.fastq.gz	/data/raw/sampleB_R2.fastq.gz	HLA-A*01:01,HLA-A*03:01,HLA-B*08:01,HLA-B*35:01,HLA-C*04:01,HLA-C*07:01	/home/faramir/repos/neondisco-jet/metadataDir/hg38.repeatmasker_reformat.txt
```

Notes:

- Use **real tabs**. In `nano`, press `Ctrl+V` then `Tab`.
- `.fastq.gz` is read directly by STAR (`--readFilesCommand zcat`); nothing is decompressed to disk. Plain `.fastq` also works. The path in the file must match what exists on disk exactly.
- HLA alleles are preferably written as `HLA-A*11:01`. The older `HLA-A11:01` form is converted automatically in Step 3.
- Alleles that pVACbind does not know (e.g. null alleles such as `HLA-B*44:309N`) are skipped with a message; the run continues with the remaining alleles.
- Lines starting with `#` are comments. A missing final newline is handled.
- Avoid sample names where one is a prefix of another followed by an underscore (e.g. `379T` and `379T_rep2`): the master script finds existing output folders with the pattern `<name>_*`.
- If you only run Step 1, the HLA and RepeatMasker columns can hold a placeholder such as `NA`.

### 6.2 HLA typing

HLA types must come from a separate tool. The aligned BAM from Step 1 (`<sample>_Aligned.sortedByCoord.out.bam` + `.bai`) is a suitable input for RNA-based typers such as arcasHLA (with HLA-HD as a fallback). A pVACbind run takes up to 6 class I alleles (2 × A, B and C).

---

## 7. Pre-Run Checklist

### 7.1 Make scripts executable

```bash
cd ~/repos/neondisco-jet/JET_identification_pipeline
chmod +x run_TEAL_pipeline.sh step1_alignment.sh step2_jet_identification.sh \
         step3_pvacbind.sh step4_compiler.sh check_integrity.sh
```

### 7.2 Check the tools

```bash
~/repos/neondisco-jet/.pixi/envs/default/bin/STAR --version
~/repos/neondisco-jet/.pixi/envs/default/bin/samtools --version | head -1
/opt/R/4.5.3/bin/Rscript --version
apptainer --version
```

### 7.3 Check the pVACtools image

```bash
SIF=/home/faramir/images/pVACtools.sif

apptainer exec $SIF pvactools --version
apptainer exec $SIF ls /opt/iedb                                   # expect: LICENSE  mhc_i  mhc_ii
apptainer exec $SIF pvacbind run --help | grep -oE "NetMHCpanEL|MHCflurryEL|BigMHC_IM|PRIME|DeepImmuno" | sort -u
```

Optional one-peptide smoke test (should finish in under a minute and write a `MHC_Class_I/` folder):

```bash
echo -e ">test\nAEADPLRAV" > /tmp/test.fasta
apptainer exec -B /home/faramir:/home/faramir -B /tmp:/tmp $SIF \
    pvacbind run /tmp/test.fasta TEST "HLA-A*02:01" \
    NetMHCpan NetMHCpanEL MHCflurryEL BigMHC_IM PRIME DeepImmuno \
    /tmp/pvacbind_test -e1 9 -t 1
```

### 7.4 Check the R packages

```bash
/opt/R/4.5.3/bin/Rscript -e "
  lib = '/home/faramir/repos/neondisco-jet/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu'
  .libPaths(lib)
  pkgs = c('getopt','ggplot2','data.table','scales','GenomicAlignments','GenomicFeatures',
           'GenomeInfoDb','txdbmaker','AnnotationDbi','biovizBase','ggbio','BSgenome.Hsapiens.UCSC.hg38')
  missing = pkgs[!pkgs %in% installed.packages(lib)[,'Package']]
  if(length(missing)) cat('MISSING:', paste(missing, collapse=', '), '\n') else cat('All packages OK\n')
"
```

### 7.5 Check `samples.tsv`

```bash
awk -F'\t' '!/^#/ && NF>0 {print NF, $1}' ~/repos/neondisco-jet/infoDir/samples.tsv   # every line should start with 5

while IFS=$'\t' read -r name r1 r2 hla rep || [ -n "$name" ]; do
    [[ -z "$name" || "$name" =~ ^# ]] && continue
    echo -n "$name: "
    [ -f "$r1" ]  && echo -n "R1 OK "  || echo -n "R1 MISSING "
    [ -f "$r2" ]  && echo -n "R2 OK "  || echo -n "R2 MISSING "
    [ -f "$rep" ] && echo "RepeatMasker OK" || echo "RepeatMasker MISSING"
done < ~/repos/neondisco-jet/infoDir/samples.tsv
```

---

## 8. Running the Pipeline

Run inside `tmux` so the job survives a dropped SSH connection:

```bash
tmux new -s teal          # detach: Ctrl+B, D      reattach: tmux attach -t teal
```

### 8.1 Full run (Steps 1–4, all samples)

```bash
cd ~/repos/neondisco-jet/JET_identification_pipeline

./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --genome hg38 \
    --genome-fasta /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf \
    --star-index /home/faramir/data/neondisco-jet/starIndexesDir \
    --star-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin \
    --samtools-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin \
    --rscript-bin /opt/R/4.5.3/bin/Rscript \
    --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R \
    --pvactools-sif /home/faramir/images/pVACtools.sif \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --tmp-dir /home/faramir/data/neondisco-jet/tmpData \
    --threads 16
```

One line:

```
cd ~/repos/neondisco-jet/JET_identification_pipeline && ./run_TEAL_pipeline.sh --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv --genome hg38 --genome-fasta /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf --star-index /home/faramir/data/neondisco-jet/starIndexesDir --star-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin --samtools-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin --rscript-bin /opt/R/4.5.3/bin/Rscript --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R --pvactools-sif /home/faramir/images/pVACtools.sif --outputs-dir /home/faramir/data/neondisco-jet/outputData --tmp-dir /home/faramir/data/neondisco-jet/tmpData --threads 16
```

### 8.2 Skipping steps

`--skip-step1` … `--skip-step4` can be combined. A skipped step still needs its **output** to exist; the master script then picks the **most recently modified** `<sample>_*` folder in `--outputs-dir`. If you keep several dated runs of the same sample, make sure the one you want is the newest. Each fresh Step 1 creates a new `<sample>_<YYYYMMDD>` folder.

| Goal | Flags to add | Flags you can leave out |
|---|---|---|
| Step 1 only | `--skip-step2 --skip-step3 --skip-step4` | `--rscript-bin --r-script --pvactools-sif` |
| Steps 2–4 | `--skip-step1` | `--genome-fasta --star-index --star-bin --samtools-bin` |
| Steps 3–4 | `--skip-step1 --skip-step2` | everything except `--samples --pvactools-sif --outputs-dir --threads` |
| Step 4 only | `--skip-step1 --skip-step2 --skip-step3` | everything except `--samples --outputs-dir` |

Example — Step 1 only:

```bash
./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --genome hg38 \
    --genome-fasta /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf \
    --star-index /home/faramir/data/neondisco-jet/starIndexesDir \
    --star-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin \
    --samtools-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --tmp-dir /home/faramir/data/neondisco-jet/tmpData \
    --threads 16 \
    --skip-step2 --skip-step3 --skip-step4
```

Example — resume at Step 2 once HLA types are in `samples.tsv`:

```bash
./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --genome hg38 \
    --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf \
    --rscript-bin /opt/R/4.5.3/bin/Rscript \
    --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R \
    --pvactools-sif /home/faramir/images/pVACtools.sif \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --threads 16 \
    --skip-step1
```

Example — rerun Step 3 only. Delete the old `pvacbind_size*` folders of that sample first, otherwise stale files remain:

```bash
rm -rf /home/faramir/data/neondisco-jet/outputData/<sample>_<date>/pvacbind_size{8,9,10,11}

./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --pvactools-sif /home/faramir/images/pVACtools.sif \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --threads 16 \
    --skip-step1 --skip-step2 --skip-step4
```

### 8.3 All flags

| Flag | Needed for | Default | Description |
|---|---|---|---|
| `--samples FILE` | all | — | Path to `samples.tsv` |
| `--genome BUILD` | all | `hg38` | `hg19`, `hg38`, `mm9` or `mm10` |
| `--genome-fasta PATH` | Step 1 | — | Genome FASTA |
| `--gtf PATH` | Steps 1–2 | — | Gene annotation GTF |
| `--star-index DIR` | Step 1 | — | STAR index directory (built if missing) |
| `--star-bin DIR` | Step 1 | — | Directory containing `STAR` |
| `--samtools-bin DIR` | Step 1 | — | Directory containing `samtools` |
| `--rscript-bin PATH` | Step 2 | — | Path to `Rscript` |
| `--r-script PATH` | Step 2 | — | Path to `JET_analysis_filtered.R` |
| `--pvactools-sif PATH` | Step 3 | — | Path to the pVACtools Apptainer image |
| `--outputs-dir DIR` | all | `./outputData` | Root output directory |
| `--tmp-dir DIR` | all | `./tmpData` | STAR temporary directory |
| `--threads N` | Steps 1, 3 | `8` | Threads for STAR and pVACbind |
| `--read-length N` | Step 1 | `100` | Read length, used as STAR `sjdbOverhang` |
| `--skip-step1` … `--skip-step4` | — | false | Skip that step |
| `-h`, `--help` | — | — | Show help |

Environment variable: `KEEP_TMP=true ./run_TEAL_pipeline.sh …` keeps pVACbind's intermediate files (default: removed).

### 8.4 Workflow for a new cohort (no HLA typing yet)

1. Create `samples.tsv` with `NA` in the HLA and RepeatMasker columns.
2. Run Step 1 only (see 8.2).
3. Run HLA typing on each `<sample>_Aligned.sortedByCoord.out.bam`.
4. Put the typed alleles (`HLA-A*xx:xx,…`) and the RepeatMasker path into `samples.tsv`.
5. Run Steps 2–4 with `--skip-step1`.

### 8.5 Rough run times (16 threads)

| Stage | Observed |
|---|---|
| STAR index build | ~40 min, once |
| Step 1, per sample | ~15–40 min, depending on library size |
| Step 3, per sample | ~45 min for about 1,000 peptides per size (sizes 8–11) |
| Step 4 | seconds |

---

## 9. Expected Outputs

Each sample gets its own dated folder:

```
outputData/
├── logs/
│   └── pipeline_<YYYYMMDD_HHMMSS>.log
└── <sample>_<YYYYMMDD>/
    ├── <sample>_Aligned.sortedByCoord.out.bam         # Step 1
    ├── <sample>_Aligned.sortedByCoord.out.bam.bai
    ├── <sample>_Chimeric.out.junction
    ├── <sample>_SJ.out.tab
    ├── <sample>_Log.final.out                         # alignment statistics
    ├── <sample>_ReadsPerGene.out.tab
    │
    ├── <sample>_Chimeric.out.annotatedJET.txt         # Step 2: annotated junction table
    ├── <sample>_Fusions_chim2.junc2.size{8,9,10,11}.fasta        # fusion peptides
    ├── <sample>_Fusions_chim2.junc2.size{8,9,10,11}.genomic.txt  # donor/acceptor nucleotides
    ├── <sample>_Fusions_chim2.junc2.size{8,9,10,11}.ids.txt      # peptide ID → junction
    ├── <sample>_Fusions_chim2.junc2.size{8,9,10,11}.RData
    │
    ├── pvacbind_size{8,9,10,11}/                      # Step 3
    │   └── MHC_Class_I/
    │       ├── <sample>.MHC_I.all_epitopes.tsv        # every peptide × allele prediction
    │       ├── <sample>.MHC_I.filtered.tsv            # after pVACbind's filters
    │       ├── <sample>.MHC_I.all_epitopes.aggregated.tsv
    │       └── log/inputs.yml                         # exact pVACbind parameters used
    │
    └── <sample>_final_combined_data_TEAL.tsv          # Step 4: final table
```

The master log records every step with timestamps. Some lines appear twice in the console and in the log; this is cosmetic.

---

## 10. Interpreting the Results

### 10.1 `annotatedJET.txt` (Step 2)

One row per chimeric or spliced junction found by STAR. Important columns:

| Column | Meaning |
|---|---|
| `file` | Where the junction came from: `Chimeric` (STAR chimeric output) or `Junction` (`SJ.out.tab`) |
| `n` | Supporting reads. Junctions are kept only if `n / uniquely mapped reads ≥ 2×10⁻⁷` |
| `Donor`, `Acceptor` | Three flags `R;E;P` per side: overlaps a **R**epeat, **E**xon (CDS) or **P**romoter (±2 kb) |
| `filter` | `TRUE` when one side is a repeat and the other side an exon (`filter.gene`) or promoter (`filter.promoter`). These junctions are the JETs that are translated. |

### 10.2 `.ids.txt` and the `Index` key

Each FASTA entry has a numeric ID that restarts at 1 for every peptide size. In `ids.txt` it is mapped to one or more entries joined with `/`:

```
457	R1R2;ORF1;chr19:44818637:+:ENST00000270233>chr19:44818713:+:(CCT)n;width=13/R1R2;ORF1;...
```

Each entry is `order;ORF;donor>acceptor;width=N`:

- `donor>acceptor` is `chr:position:strand:feature` for both ends. The feature is a transcript ID for the exon side and a repeat name for the TE side. One peptide can map to several transcripts of the same gene.
- `width=N` is the number of amino acids of the translated reading frame (`ORF1`–`ORF3`) after truncation at the first stop codon. It is smaller than the translated window when a stop codon ends the frame inside the junction region.

pVACbind's `Index` column equals this FASTA ID, so **`(Size, Index)` is the key** that links a prediction to its junction.

### 10.3 pVACbind outputs (Step 3)

- `all_epitopes.tsv` — every sub-peptide × HLA allele, all algorithm scores. Nothing filtered.
- `filtered.tsv` — the same columns after pVACbind's binding filter and top-score filter (both steps appear in the pVACbind log). Check the exact thresholds in `MHC_Class_I/log/inputs.yml`; in the runs so far they were IC50 < 500 nM, percentile < 2 and `top_score_metric: median`.
- The top-score filter normally reduces the table to the best-scoring epitope per FASTA entry. If you need every binding window and allele, compile from `all_epitopes.tsv` instead. A quick check: compare `grep -c '^>' <fasta>` with `tail -n +2 <filtered.tsv> | wc -l`; if the table has at most one row per entry, the top-score filter is active.

Common cut-offs: IC50 < 50 nM strong, < 500 nM weak/intermediate binder; percentile rank ≤ 0.5 % strong, ≤ 2 % weak. Lower is better for IC50 and percentile ranks. For the immunogenicity predictors (BigMHC_IM, PRIME, DeepImmuno), check each tool's documentation for the direction of the score. The pVACbind columns also contain manufacturability flags (cysteine counts, N-terminal asparagine, etc.), and NetMHCstabpan stability columns when `--netmhc-stab` is on.

### 10.4 `<sample>_final_combined_data_TEAL.tsv` (Step 4)

All columns of the four `filtered.tsv` files, stacked, with three added columns:

| Position | Column | Content |
|---|---|---|
| 1 | `Sample` | Sample name |
| 2 | `Size` | Peptide length (8–11) that the row came from |
| 3 … N−1 | pVACbind columns, unchanged | `Index`, `HLA Allele`, `Sub-peptide Position`, `Epitope Seq`, `Median IC50 Score`, `Best IC50 Score`, `Best IC50 Score Method`, `Median Percentile`, `Best Percentile`, per-algorithm scores and percentiles, manufacturability flags, … |
| last | `Junction` | `donor>acceptor[width=N]`, several entries separated by `;`; `NA` if the `Index` is not found in `ids.txt` |

The exact set of pVACbind columns depends on the pVACtools version and flags.

Using it downstream:

```r
library(readr); library(dplyr)

df <- read_tsv("379T_final_combined_data_TEAL.tsv")
df %>% filter(`Best IC50 Score` < 50) %>% count(Size)
```

Combining several samples (assumes identical settings, hence identical columns):

```bash
files=(/home/faramir/data/neondisco-jet/outputData/*/*_final_combined_data_TEAL.tsv)
head -1 "${files[0]}" > all_samples_TEAL.tsv
for f in "${files[@]}"; do tail -n +2 "$f" >> all_samples_TEAL.tsv; done
```

---

## 11. Integrity Check

`check_integrity.sh` compares two runs of the same sample by md5sum and line count (BAMs by size, since their headers carry timestamps):

```bash
./check_integrity.sh -o 379T_20260701 -n 379T_20260715                  # sample name derived from -o
./check_integrity.sh -o 379T_20260701 -n 379T_20260715 -s 379T
./check_integrity.sh -o 379T_20260701 -n 379T_20260715 -d /mnt/data/neondisco-jet/outputData
```

Reports are saved to `outputData/integrity_check/`. The Step 1–2 file groups still apply. The Step 3–4 entries in `FILE_GROUPS` refer to the old netMHCpan files (`*.netmhcpan4.txt`, `all_sizes_binders.tsv`); edit them to `pvacbind_size*/MHC_Class_I/*.tsv` and `*_final_combined_data_TEAL.tsv` before using the script on pVACbind runs.

---

## 12. Moving to Another Machine

Some paths are set inside the scripts rather than on the command line:

| Where | What |
|---|---|
| `JET_analysis_filtered.R` | `lib.local` (renv library path) and `txdb.cache` (cached TxDb location); `--gtf` falls back to a default path if not given |
| `step3_pvacbind.sh` | The `apptainer exec … -B /GB:/GB -B /home/faramir:/home/faramir` bind mounts. Every directory holding inputs or outputs must be covered by a bind, or the container cannot see it. Add or change `-B` entries for your storage layout. |
| `step3_pvacbind.sh` | `--iedb-install-directory /opt/iedb` is a path **inside** the image |

---

## 13. Troubleshooting

**STAR: "genomeParameters.txt missing" or index errors**
The index build may have been interrupted. Delete the index directory and rerun; Step 1 rebuilds it. Check that `--genome-fasta` and `--gtf` are correct.

**STAR: "quality string length is not equal to sequence length"**
The FASTQ is truncated or damaged. Test the file with `gzip -t file.fastq.gz`, and check lengths with
`zcat file.fastq.gz | awk 'NR%4==2{s=length($0)} NR%4==0 && length($0)!=s{print NR; exit}'`.

**"No output directory found for sample X"**
You used `--skip-step1` but there is no `<X>_*` folder in `--outputs-dir`. Run Step 1 first.

**Step 2: "package not found"**
Check `lib.local` at the top of the R script with `Rscript -e ".libPaths()"`.

**Step 3: pVACbind prints `Allele … not valid. Skipping.`**
Expected for null or unsupported alleles. The run continues with the rest. If every allele is invalid, no output is written and Step 3 stops with an error.

**Step 3: `Epitope Length 8 is not valid for Method DeepImmuno`**
DeepImmuno supports fewer alleles and lengths than the other predictors; those combinations are skipped.

**Step 3: `NetMHCstabpan substituted HLA-… for HLA-…`**
Informational. NetMHCstabpan uses the closest allele it supports; the distance is printed.

**Step 3: thousands of `did not find allele distance info from netmhcpan_output` warnings and stray `tmp…` files in the working directory**
Caused by running the container with `--no-mount tmp`. Do not add that option; leave `/tmp` mounted.

**Step 3: very slow, log shows `IEDB API Output doesn't match input. Retrying.`**
pVACbind fell back to the IEDB web API. Make sure the command contains `--iedb-install-directory /opt/iedb` and that `ls /opt/iedb` works inside the image.

**Step 3: `UserWarning: pkg_resources is deprecated`**
Harmless; it comes from the IEDB code inside the image.

**Step 3: "No such file" for an input inside the container**
The file's directory is not bound. Add a `-B` entry (Section 12).

**Step 4: "No pVACbind filtered.tsv + ids.txt found"**
Step 3 has not completed for that sample. A size with an empty `filtered.tsv` is skipped with a message.

---

## 14. Citation

If you use this pipeline, please cite:

> Burbage M, Rocañín-Arjó A, et al. Epigenetically-controlled tumor antigens derived from splice junctions between exons and transposable elements. **Science Immunology**, 2023.

> Hundal J, et al. pVACtools: a computational toolkit to identify and visualize cancer neoantigens. **Cancer Immunology Research**, 2020.

The individual predictors run through pVACbind (NetMHCpan, MHCflurry, BigMHC, PRIME, DeepImmuno, NetMHCstabpan) have their own publications and licences; see the pVACtools documentation.

---

## 15. Authors and Contact

**Original pipeline:** Alexandre Houy and Christel Goudot (Institut Curie); modified by Ares Rocañín-Arjó (07/2022).

**Adapted for human GRCh38 and extended (pVACbind, multi-sample orchestration) by:** Muhammad Faramir — Cancer Research Malaysia (CRMY).

**Contact for the original method:**
- Ares Rocañín-Arjó: maria-ares.rocanin-arjo@curie.fr
- Marianne Burbage: marianne.burbage@curie.fr
- Christel Goudot: christel.goudot@curie.fr
