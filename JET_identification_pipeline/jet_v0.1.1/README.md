# NeonDisco-TEAL — JET Identification Pipeline

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

> **JET Identification Pipeline — adapted for human GRCh38**
> Burbage & Rocañín-Arjó et al., *Science Immunology* 2023

---

## Table of Contents

1. [Background](#1-background)
2. [Pipeline Overview](#2-pipeline-overview)
3. [Repository Structure](#3-repository-structure)
4. [Software Requirements](#4-software-requirements)
5. [Reference Files Required](#5-reference-files-required)
6. [Input Data Preparation](#6-input-data-preparation)
7. [Directory Setup](#7-directory-setup)
8. [Pre-Run Preparation](#8-pre-run-preparation)
9. [Running the Pipeline](#9-running-the-pipeline)
10. [Expected Outputs](#10-expected-outputs)
11. [Interpreting the Results](#11-interpreting-the-results)
12. [Integrity Check](#12-integrity-check)
13. [Troubleshooting](#13-troubleshooting)
14. [Citation](#14-citation)
15. [Authors and Contact](#15-authors-and-contact)

---

## 1. Background

Transposable elements (TEs) make up nearly half of the human genome and are normally silenced by epigenetic mechanisms such as H3K9me3 deposition via SETDB1. In tumour cells, this epigenetic silencing is frequently disrupted, allowing TEs to become aberrantly transcribed and spliced into gene transcripts.

**JETs (Junctions between Exons and Transposable Elements)** are non-canonical splice junctions that arise when an exonic donor or acceptor splice site is paired with a cryptic splice site within an intronic or intergenic TE. The resulting fusion transcript encodes a chimeric peptide — part gene product, part transposable element — that is entirely absent from normal tissues. Because these peptides are presented on MHC class I molecules at the tumour cell surface, they represent potential tumour-specific neoantigens that could be recognised by cytotoxic T cells.

This pipeline identifies, quantifies, and selects JETs from RNA-seq data and predicts which resulting fusion peptides are likely to bind MHC class I molecules. It was originally developed by Alexandre Houy and Christel Goudot at Institut Curie and described in:

> Burbage M\*, Rocañín-Arjó A\* et al. *Epigenetically controlled tumor antigens derived from splice junctions between exons and transposable elements.* **Science Immunology**, 2023.

This repository adapts the original pipeline for **human GRCh38** samples and wraps all steps into a single master orchestrator (`run_TEAL_pipeline.sh`) that supports multiple samples in a single run.

---

## 2. Pipeline Overview

The pipeline consists of four sequential steps:

```
RNA-seq FASTQ
      │
      ▼
┌─────────────────────────────────────────────┐
│  STEP 1 — STAR 2-pass Alignment             │
│  Auto-builds genome index if not found;     │
│  identifies chimeric reads spanning         │
│  exon-TE boundaries; sorts and indexes BAMs │
└────────────────────┬────────────────────────┘
                     │  Chimeric.out.junction
                     │  SJ.out.tab
                     ▼
┌─────────────────────────────────────────────┐
│  STEP 2 — JET Identification (R script)     │
│  Annotates junctions against RepeatMasker   │
│  and gene annotations; filters for JETs;    │
│  extracts fusion peptide sequences          │
└────────────────────┬────────────────────────┘
                     │  .fasta (sizes 8–11)
                     ▼
┌─────────────────────────────────────────────┐
│  STEP 3 — MHC-I Binding Prediction          │
│  netMHCpan 4.x predicts binding affinity    │
│  for each fusion peptide against patient    │
│  HLA alleles                                │
└────────────────────┬────────────────────────┘
                     │  .netmhcpan4.txt (sizes 8–11)
                     ▼
┌─────────────────────────────────────────────┐
│  STEP 4 — Binder Aggregation                │
│  Merges results across peptide sizes 8–11;  │
│  reports binders (rank<2%) and strong       │
│  binders (rank<0.5%) with junction mapping  │
└────────────────────┬────────────────────────┘
                     │
                     ▼
              all_sizes_binders.tsv
```

---

## 3. Repository Structure

```
neondisco-jet/
├── JET_identification_pipeline/
│   ├── run_TEAL_pipeline.sh              # Master orchestrator (entry point)
│   ├── step1_alignment.sh               # STAR alignment + BAM processing
│   ├── step2_jet_identification.sh      # JET calling wrapper (calls R script)
│   ├── step3_netmhcpan.sh              # netMHCpan binding prediction
│   ├── aggregate_binders.sh            # Cross-size binder aggregation
│   ├── JET_analysis_filtered.R         # Core R script for JET annotation
│   ├── check_integrity.sh              # md5sum integrity checker
│   └── lib/
│       └── common.sh                   # Shared logging functions
├── infoDir/
│   ├── samples.txt                     # Legacy 3-column metadata (original scripts)
│   └── samples.tsv                     # New 5-column metadata (pipeline)
├── inputData/
│   ├── <sample>_R1.fastq(.gz)         # RNA-seq reads
│   ├── <sample>_R2.fastq(.gz)
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   └── Homo_sapiens.GRCh38.115.gtf
├── metadataDir/
│   ├── <sample>.repeatmasker_reformat.txt
│   └── GRCh38_115.txdb                # Cached TxDb (auto-generated on first run)
├── dependencies/
│   └── netMHCpan-4.2/                 # netMHCpan installation
└── renv/                              # R package library (managed by renv)
```

---

## 4. Software Requirements

| Tool | Version | How to install |
|---|---|---|
| STAR | 2.5.3a | `pixi add star=2.5.3a` |
| samtools | 1.3 | `pixi add samtools=1.3` |
| R | 4.5.x | via `rig` (system-wide) |
| netMHCpan | 4.1 or 4.2 | Download from DTU Health Tech (license required) |

### R packages required

All packages are managed via `renv` and stored in the project library. Install from within the project directory (`~/repos/neondisco-jet/`):

```r
# In R, from ~/repos/neondisco-jet/
BiocManager::install(c(
    "GenomicAlignments",
    "GenomicFeatures",
    "GenomeInfoDb",
    "biovizBase",
    "ggbio",
    "BSgenome.Hsapiens.UCSC.hg38",
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
))
install.packages(c("ggplot2", "data.table", "scales", "getopt"))
```

> **Note:** `biovizBase` requires the system library `libuv1-dev`. If installation fails, ask your sysadmin to run `sudo apt install libuv1-dev`.

### pixi environment setup (STAR + samtools)

```bash
cd ~/repos/neondisco-jet
pixi add star=2.5.3a samtools=1.3
```

### netMHCpan installation

netMHCpan requires a free academic license from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/). After downloading:

```bash
mkdir -p ~/repos/neondisco-jet/dependencies
tar -xzf netMHCpan-4.2.Linux.tar.gz -C ~/repos/neondisco-jet/dependencies/
# Follow the configuration instructions in the netMHCpan README
```

---

## 5. Reference Files Required

### 5.1 Genome FASTA

Download the GRCh38 primary assembly from Ensembl:

```bash
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/repos/neondisco-jet/inputData/
```

### 5.2 Gene Annotation GTF

```bash
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
mv Homo_sapiens.GRCh38.115.gtf ~/repos/neondisco-jet/inputData/
```

### 5.3 RepeatMasker Annotation (per sample, from UCSC Table Browser)

This is a **sample-specific** file that must be downloaded from UCSC for the relevant genome. It is used in Step 2 to annotate which junctions overlap transposable elements.

**Steps to download:**

1. Go to [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
2. Set the following options:

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
| Output file | `<sample>.repeatmasker_reformat.txt` |
| File type returned | plain text |

3. Click **Get output**, then select these columns only:
   - `genoName`
   - `genoStart`
   - `genoEnd`
   - `strand`
   - `repName`
   - `repClass`
   - `repFamily`

4. Click **Get output** and save the file.

5. Move it to:
```bash
mv <sample>.repeatmasker_reformat.txt ~/repos/neondisco-jet/metadataDir/
```

> **Note:** The RepeatMasker file uses UCSC chromosome naming (`chr1`, `chr2`, ...). The pipeline handles the conversion from Ensembl naming (`1`, `2`, ...) automatically in the R script.

### 5.4 STAR Genome Index

**You do not need to build the index manually.** `step1_alignment.sh` checks for the index automatically at the start of every run:

- If the index **does not exist** → it builds it automatically before alignment. The index directory is also created automatically via `mkdir -p`. This takes approximately 30–40 minutes and requires ~30 GB of RAM on the first run only.
- If the index **already exists** (detected by the presence of `genomeParameters.txt` inside `--star-index`) → it skips the build entirely and proceeds straight to alignment.

All you need to do is pass the desired index directory path via `--star-index` when running the pipeline. The same path is reused for all subsequent samples — the index only ever builds once per genome.

```bash
# Example: first run on a fresh server — index will be built automatically
./run_TEAL_pipeline.sh \
    --star-index /home/faramir/data/neondisco-jet/starIndexesDir \
    ...

# Subsequent runs — index already exists, build is skipped automatically
./run_TEAL_pipeline.sh \
    --star-index /home/faramir/data/neondisco-jet/starIndexesDir \
    ...
```

> **Note:** The index is genome-level and shared across all samples. On a shared server, coordinate with other users before triggering a rebuild to avoid duplicate index builds consuming disk and RAM simultaneously.

---

## 6. Input Data Preparation

### 6.1 Samples metadata file (`samples.tsv`)

Create a tab-separated file with one row per sample and **5 columns** in this exact order:

```
name<TAB>R1_path<TAB>R2_path<TAB>HLA_alleles<TAB>repeatmasker_file
```

| Column | Description | Example |
|---|---|---|
| `name` | Unique sample identifier | `379T` |
| `R1_path` | Absolute path to R1 FASTQ | `/home/faramir/repos/neondisco-jet/inputData/379T_R1.fastq` |
| `R2_path` | Absolute path to R2 FASTQ | `/home/faramir/repos/neondisco-jet/inputData/379T_R2.fastq` |
| `HLA_alleles` | Comma-separated HLA alleles in netMHCpan format | `HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04` |
| `repeatmasker_file` | Absolute path to sample RepeatMasker file | `/home/faramir/repos/neondisco-jet/metadataDir/379T.repeatmasker_reformat.txt` |

**Example `samples.tsv`:**

```
379T	/home/faramir/repos/neondisco-jet/inputData/379T_R1.fastq	/home/faramir/repos/neondisco-jet/inputData/379T_R2.fastq	HLA-A11:01,HLA-A02:03,HLA-B13:01,HLA-B38:02,HLA-C07:02,HLA-C03:04	/home/faramir/repos/neondisco-jet/metadataDir/379T.repeatmasker_reformat.txt
319T	/home/faramir/data/raw/319T_R1.fastq	/home/faramir/data/raw/319T_R2.fastq	HLA-A02:01,HLA-A24:02,HLA-B40:01,HLA-B51:01,HLA-C03:03,HLA-C15:02	/home/faramir/repos/neondisco-jet/metadataDir/319T.repeatmasker_reformat.txt
```

> **Important notes:**
> - Use **real tabs** between columns (not spaces). In `nano`, press `Ctrl+V` then `Tab` to insert a literal tab.
> - FASTQ files can be either `.fastq` or `.fastq.gz`. The pipeline will auto-decompress `.gz` files if needed.
> - Lines starting with `#` are treated as comments and skipped.
> - HLA alleles must be in standard netMHCpan format: `HLA-A11:01` (not `HLA-A*11:01`).
> - Each sample name must be unique as it becomes the output folder prefix.

Save the file to:
```bash
~/repos/neondisco-jet/infoDir/samples.tsv
```

### 6.2 HLA allele typing

HLA alleles must be determined separately (e.g. via arcasHLA or OptiType on WES/WGS data). The pipeline accepts up to 6 alleles (2× HLA-A, 2× HLA-B, 2× HLA-C), which is the typical diploid set for class I.

---

## 7. Directory Setup

Create all required directories before the first run:

```bash
# Output and working directories
mkdir -p ~/data/neondisco-jet/outputData
mkdir -p ~/data/neondisco-jet/outputData/logs
mkdir -p ~/data/neondisco-jet/tmpData
mkdir -p ~/data/neondisco-jet/starIndexesDir

# Repo metadata and input directories
mkdir -p ~/repos/neondisco-jet/metadataDir
mkdir -p ~/repos/neondisco-jet/inputData
mkdir -p ~/repos/neondisco-jet/infoDir
mkdir -p ~/repos/neondisco-jet/dependencies
```

> **Directories created automatically by the pipeline (do not need to be created manually):**
>
> | Directory | Created by | When |
> |---|---|---|
> | `starIndexesDir/` | `step1_alignment.sh` | Before genome indexing on first run |
> | `outputData/<sample>_<date>/` | `step1_alignment.sh` | At the start of each sample run |
> | `outputData/logs/` | `run_TEAL_pipeline.sh` | At pipeline startup |
> | `tmpData/STAR_<sample>_<date>/` | `step1_alignment.sh` | During STAR alignment |

---

## 8. Pre-Run Preparation

### 8.1 Make all scripts executable

```bash
cd ~/repos/neondisco-jet/JET_identification_pipeline

chmod +x run_TEAL_pipeline.sh
chmod +x step1_alignment.sh
chmod +x step2_jet_identification.sh
chmod +x step3_netmhcpan.sh
chmod +x aggregate_binders.sh
chmod +x check_integrity.sh
```

### 8.2 Verify tool availability

```bash
# Check STAR
~/repos/neondisco-jet/.pixi/envs/default/bin/STAR --version
# Expected: STAR_2.5.3a

# Check samtools
~/repos/neondisco-jet/.pixi/envs/default/bin/samtools --version
# Expected: samtools 1.3

# Check Rscript
/opt/R/4.5.3/bin/Rscript --version

# Check netMHCpan
~/repos/neondisco-jet/dependencies/netMHCpan-4.2/Linux_x86_64/bin/netMHCpan -h
```

### 8.3 Verify R packages are installed

```bash
/opt/R/4.5.3/bin/Rscript -e "
  lib = '/home/faramir/repos/neondisco-jet/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu'
  .libPaths(lib)
  pkgs = c('getopt','ggplot2','data.table','scales','GenomicAlignments',
           'GenomicFeatures','GenomeInfoDb','biovizBase','ggbio',
           'BSgenome.Hsapiens.UCSC.hg38','TxDb.Hsapiens.UCSC.hg38.knownGene')
  missing = pkgs[!pkgs %in% installed.packages(lib)[,'Package']]
  if(length(missing)) cat('MISSING:', paste(missing, collapse=', '), '\n') else cat('All packages OK\n')
"
```

### 8.4 Verify your samples.tsv

```bash
# Check it is tab-separated and has 5 columns
awk -F'\t' '{print NF, $1}' ~/repos/neondisco-jet/infoDir/samples.tsv
# Every non-comment line should print "5 <sample_name>"

# Check all input files exist
while IFS=$'\t' read -r name r1 r2 hla rep; do
    [[ "$name" =~ ^# ]] && continue
    echo -n "Sample ${name}: "
    [ -f "$r1" ] && echo -n "R1 OK " || echo -n "R1 MISSING "
    [ -f "$r2" ] && echo -n "R2 OK " || echo -n "R2 MISSING "
    [ -f "$rep" ] && echo "RepMask OK" || echo "RepMask MISSING"
done < ~/repos/neondisco-jet/infoDir/samples.tsv
```

---

## 9. Running the Pipeline

### 9.1 Full run (all 4 steps)
Multiple lines:

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
    --netmhcpan-bin /home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2 \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --tmp-dir /home/faramir/data/neondisco-jet/tmpData \
    --threads 16
```

One liner:

```
cd ~/repos/neondisco-jet/JET_identification_pipeline && ./run_TEAL_pipeline.sh --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv --genome hg38 --genome-fasta /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf --star-index /home/faramir/data/neondisco-jet/starIndexesDir --star-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin --samtools-bin /home/faramir/repos/neondisco-jet/.pixi/envs/default/bin --rscript-bin /opt/R/4.5.3/bin/Rscript --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R --netmhcpan-bin /home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2 --outputs-dir /home/faramir/data/neondisco-jet/outputData --tmp-dir /home/faramir/data/neondisco-jet/tmpData --threads 16
```

It is strongly recommended to run this inside a `tmux` session so it continues if your SSH connection drops:

```bash
tmux new -s jet_pipeline
# then run the command above
# detach with Ctrl+B, D
# reattach later with: tmux attach -t jet_pipeline
```

### 9.2 Skip steps (for reprocessing)

If Step 1 (STAR alignment) has already been run successfully and you only want to rerun Steps 2–4:
Multiple lines:
```bash
./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --genome hg38 \
    --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf \
    --rscript-bin /opt/R/4.5.3/bin/Rscript \
    --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R \
    --netmhcpan-bin /home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2/Linux_x86_64/bin \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --skip-step1
```
One liner:

```
cd ~/repos/neondisco-jet/JET_identification_pipeline && ./run_TEAL_pipeline.sh --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv --genome hg38 --gtf /home/faramir/repos/neondisco-jet/inputData/Homo_sapiens.GRCh38.115.gtf --rscript-bin /opt/R/4.5.3/bin/Rscript --r-script /home/faramir/repos/neondisco-jet/JET_identification_pipeline/JET_analysis_filtered.R --netmhcpan-bin /home/faramir/repos/neondisco-jet/dependencies/netMHCpan-4.2 --outputs-dir /home/faramir/data/neondisco-jet/outputData --skip-step1
```

To rerun only the binder aggregation (e.g. after changing thresholds):
Multiple lines:

```bash
./run_TEAL_pipeline.sh \
    --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv \
    --outputs-dir /home/faramir/data/neondisco-jet/outputData \
    --skip-step1 --skip-step2 --skip-step3
```

One liner:

```
cd ~/repos/neondisco-jet/JET_identification_pipeline && ./run_TEAL_pipeline.sh --samples /home/faramir/repos/neondisco-jet/infoDir/samples.tsv --outputs-dir /home/faramir/data/neondisco-jet/outputData --skip-step1 --skip-step2 --skip-step3
```

### 9.3 All available flags

| Flag | Required | Default | Description |
|---|---|---|---|
| `--samples FILE` | ✅ | — | Path to samples.tsv |
| `--genome BUILD` | | `hg38` | Genome build: `hg19`, `hg38`, `mm9`, `mm10` |
| `--genome-fasta PATH` | ✅ (Step 1) | — | Genome FASTA file |
| `--gtf PATH` | ✅ (Steps 1–2) | — | Gene annotation GTF |
| `--star-index DIR` | ✅ (Step 1) | — | STAR genome index directory — built automatically on first run if missing |
| `--star-bin DIR` | ✅ (Step 1) | — | Directory containing `STAR` binary |
| `--samtools-bin DIR` | ✅ (Step 1) | — | Directory containing `samtools` binary |
| `--rscript-bin PATH` | ✅ (Step 2) | — | Path to `Rscript` binary |
| `--r-script PATH` | ✅ (Step 2) | — | Path to `JET_analysis_filtered.R` |
| `--netmhcpan-bin DIR` | ✅ (Step 3) | — | Directory containing `netMHCpan` binary |
| `--outputs-dir DIR` | | `./outputData` | Root output directory |
| `--tmp-dir DIR` | | `./tmpData` | STAR temporary directory |
| `--threads N` | | `8` | Number of threads for STAR/samtools |
| `--read-length N` | | `100` | RNA-seq read length (used as `sjdbOverhang`) |
| `--skip-step1` | | false | Skip STAR alignment |
| `--skip-step2` | | false | Skip JET identification |
| `--skip-step3` | | false | Skip netMHCpan |
| `--skip-step4` | | false | Skip binder aggregation |
| `-h`, `--help` | | — | Show help message |

---

## 10. Expected Outputs

For each sample, all outputs are written to a timestamped directory:

```
outputData/
└── <sample>_<YYYYMMDD>/
    ├── <sample>_Aligned.sortedByCoord.out.bam          # Genome-aligned BAM
    ├── <sample>_Aligned.sortedByCoord.out.bam.bai      # BAM index
    ├── <sample>_Chimeric.out.junction                  # Raw STAR chimeric junctions
    ├── <sample>_Chimeric.out.sort.bam                  # Chimeric reads BAM
    ├── <sample>_Chimeric.out.sort.bam.bai              # Chimeric BAM index
    ├── <sample>_SJ.out.tab                             # Splice junction table
    ├── <sample>_Log.final.out                          # STAR alignment statistics
    ├── <sample>_ReadsPerGene.out.tab                   # Gene-level read counts
    │
    ├── <sample>_Chimeric.out.annotatedJET.txt          # Annotated junction table (Step 2)
    │
    ├── <sample>_Fusions_chim2.junc2.size8.fasta        # Fusion peptide sequences (size 8)
    ├── <sample>_Fusions_chim2.junc2.size8.genomic.txt  # Donor/acceptor nucleotide sequences
    ├── <sample>_Fusions_chim2.junc2.size8.ids.txt      # Peptide ID → junction mapping
    ├── <sample>_Fusions_chim2.junc2.size8.netmhcpan4.txt  # netMHCpan output
    ├── <sample>_Fusions_chim2.junc2.size8.RData        # Full R workspace
    │   ... (same files repeated for sizes 9, 10, 11)
    │
    └── all_sizes_binders.tsv                           # Final aggregated binder table (Step 4)
```

A pipeline log is written to:

```
outputData/
└── logs/
    └── pipeline_<YYYYMMDD_HHMMSS>.log
```

---

## 11. Interpreting the Results

### 11.1 `annotatedJET.txt`

The main junction annotation table produced by Step 2. Each row is one chimeric junction with columns indicating whether the donor and acceptor sides overlap a repeat element (R), exon (E), or promoter (P). Rows where one side is `R` and the other is `E` are the true JETs.

### 11.2 `all_sizes_binders.tsv`

The final output of the full pipeline. Columns:

| Column | Description |
|---|---|
| `Size` | Peptide length tested (8, 9, 10, or 11 amino acids) |
| `ID` | Internal FASTA sequence ID (cross-reference with `ids.txt` for that size) |
| `Peptide` | Amino acid sequence of the predicted binder |
| `MinRank` | Best (lowest) percentile rank across all HLA alleles — the key filter column |
| `NB` | Number of HLA alleles predicted to bind this peptide (netMHCpan internal count) |
| `Junction` | Donor>Acceptor junction coordinates identifying the source JET |

**Filtering thresholds:**

| Threshold | Category | Recommended use |
|---|---|---|
| `MinRank < 2.0` | Binder (weak + strong) | Primary filter — used in the original paper |
| `MinRank < 0.5` | Strong binder | Stricter filter for prioritisation |

**Expected numbers (from the paper, mouse melanoma cell lines):**

- ~400–450 JETs detected per sample
- ~50% have at least one MHC-I binding peptide (rank < 2%)
- ~44–56% of JETs validated by RT-PCR in the paper

### 11.3 Quick summary from the command line

```bash
# Print binder counts for a specific sample
sample_dir=/home/faramir/data/neondisco-jet/outputData/379T_20260615

total=$(tail -n +2 "${sample_dir}/all_sizes_binders.tsv" | wc -l)
strong=$(tail -n +2 "${sample_dir}/all_sizes_binders.tsv" | awk -F'\t' '$4<0.5' | wc -l)
uniq=$(tail -n +2 "${sample_dir}/all_sizes_binders.tsv" | awk -F'\t' '{print $6}' | tr ';' '\n' | sort -u | wc -l)

echo "Total binders (rank<2)  : ${total}"
echo "Strong binders (rank<0.5): ${strong}"
echo "Unique JET junctions     : ${uniq}"
```

---

## 12. Integrity Check

To verify that a new run produces results consistent with a reference run, use the provided integrity check script:

```bash
cd ~/repos/neondisco-jet/JET_identification_pipeline

# Basic — sample name auto-derived from directory name
./check_integrity.sh -o 379T_20260701 -n 379T_20260715

# Explicit sample name
./check_integrity.sh -o 379T_20260701 -n 379T_20260715 -s 379T

# Different outputs base directory
./check_integrity.sh -o 379T_20260701 -n 379T_20260715 -d /mnt/data/neondisco-jet/outputData

# Show help
./check_integrity.sh -h
```

This script:
- Computes md5sums for all key intermediate and final files
- Compares line counts between the reference and new run
- Reports a per-file `PASS` / `WARN` / `FAIL` status
- Prints a summary of binder counts from both runs side by side
- Saves a full report to `outputData/integrity_check/`

> **Note:** BAM files are compared by file size rather than md5sum, because BAM headers embed timestamps that differ between runs even when the alignments are identical.

---

## 13. Troubleshooting

**STAR alignment fails with "genome not found" or "genomeParameters.txt missing"**
→ The index build may have failed or been interrupted. Delete the `starIndexesDir` folder and re-run the pipeline — `step1_alignment.sh` will detect the missing `genomeParameters.txt` and rebuild the index automatically. Make sure `--genome-fasta` and `--gtf` paths are correct and the files are not corrupted.

**STAR index is rebuilding on every run**
→ Check that `--star-index` points to the same directory each run and that `genomeParameters.txt` is present inside it after a successful build. If the path differs between runs the check will always trigger a rebuild.

**R script fails: `error: unable to find an inherited method for function 'seqinfo' for signature 'x = "NULL"'`**
→ The `ideo[["hg38"]]` lookup failed. The R script should use `ideo[["hg19"]]` for the ideogram object — check lines around `genome.ideo` in `JET_analysis_filtered.R`.

**R script fails: package not found**
→ Verify that `lib.local` at the top of `JET_analysis_filtered.R` points to your actual renv library path. Run `Rscript -e ".libPaths()"` to check.

**netMHCpan fails: "command not found"**
→ Confirm the binary path with `find ~/repos/neondisco-jet/dependencies -name "netMHCpan"`.

**"No output directory found for sample X"**
→ When using `--skip-step1`, the pipeline looks for an existing `outputData/<sample>_*` directory. If none exists, Step 1 must be run first.

**`samples.tsv` validation fails: "R1 fastq not found"**
→ Check that paths in `samples.tsv` are absolute and that the files exist exactly as written (`.fastq` vs `.fastq.gz` must match what is on disk).

**BAM index missing**
→ Re-run `samtools index` manually on the BAM file, or re-run Step 1 with `--skip-step2 --skip-step3 --skip-step4`.

---

## 14. Citation

If you use this pipeline, please cite:

> Burbage M\*, Rocañín-Arjó A\*, Baudon B, Goudot C, Indirect recognition of TE-derived antigens by cytotoxic T lymphocytes involves tumor immunosurveillance. **Science Immunology**, 2023.

---

## 15. Authors and Contact

**Original pipeline:**
- Alexandre Houy — Institut Curie
- Christel Goudot — Institut Curie
- Ares Rocañín-Arjó — Institut Curie

**Adapted for human GRCh38 by:**
- Muhammad Faramir — Cancer Research Malaysia (CRMY)

**Original authors contact:**
- Ares Rocañín-Arjó: maria-ares.rocanin-arjo@curie.fr
- Marianne Burbage: marianne.burbage@curie.fr
- Christel Goudot: christel.goudot@curie.fr
