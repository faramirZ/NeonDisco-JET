# NeonDisco-TEAL v.2.1.0 | JET Identification Pipeline

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
3. [Choose a Workflow](#3-choose-a-workflow)
4. [Before You Run](#4-before-you-run)
5. [Run with Nextflow](#5-run-with-nextflow)
6. [Run with Bash](#6-run-with-bash)
7. [Outputs and Run Reports](#7-outputs-and-run-reports)
8. [Citation and Credits](#8-citation-and-credits)

## 1. Background

Transposable elements (TEs) make up a large part of the human genome and are usually kept silent by epigenetic mechanisms. In some tumours, this silencing is disrupted. A transcribed TE can then be spliced together with a gene exon.

**JETs (Junctions between Exons and Transposable Elements)** are splice junctions that connect a gene exon to a transposable element. A transcript containing a JET can encode a peptide that combines sequence from the gene and the TE. If the peptide is presented by a tumour cell on an HLA class I molecule, it may be a candidate tumour antigen.

NeonDisco-TEAL adapts the original Institut Curie JET workflow for human GRCh38. It uses paired-end RNA-seq to identify JETs, generates candidate peptides of 8–11 amino acids, and uses pVACbind to predict peptide–HLA interactions for each sample. The results help prioritise candidates for further investigation.

## 2. Pipeline Overview

| Step | What it does | Main outputs |
|---|---|---|
| **1. STAR alignment** | Uses an existing STAR index or builds one from the reference, then aligns paired-end RNA-seq reads. | Sorted BAM and BAM index, chimeric junctions, splice junctions and the STAR alignment log |
| **2. JET identification** | Uses an R script, gene annotation and RepeatMasker data to identify JETs and generate peptides of lengths 8–11. | Annotated JET table, peptide FASTAs and peptide ID tables |
| **3. pVACbind** | Runs MHC class I predictions with each sample's HLA alleles. | Prediction reports by sample and peptide length |
| **4. Compile results** | Combines filtered pVACbind results and links peptides back to their junctions. This step is available in the Bash workflow. | `<sample>_final_combined_data_TEAL.tsv` |

## 3. Choose a Workflow

| Workflow | Directory and detailed guide | What it runs |
|---|---|---|
| **Nextflow** | [`NF_JET_identificatio_pipeline/`](NF_JET_identificatio_pipeline/README.md) | Steps 1–3. Runs samples and peptide lengths in parallel as resources allow. |
| **Bash** | [`JET_identification_pipeline/`](JET_identification_pipeline/README.md) | Steps 1–4. The master script processes the sample list and writes the final combined TEAL table for each sample. |

Both workflows accept multiple samples in one run. With **Nextflow**, each sample can move to the next step as soon as its own inputs are ready. For example, `379T` can start JET identification while `378T` is still aligning. Nextflow also supports running individual stages and reusing eligible completed tasks with `-resume`.

Use the detailed README inside the corresponding directory for full setup instructions, input formats, step-by-step commands and troubleshooting.

## 4. Before You Run

**Launch both workflows from the repository root**, where `pixi.toml`, `NF_JET_identificatio_pipeline/` and `JET_identification_pipeline/` are located. On the development server this is `/home/faramir/repos/neondisco-jet`. On another machine, use the folder where you cloned or extracted the repository.

```bash
git clone https://github.com/faramirZ/NeonDisco-JET.git neondisco-jet
cd neondisco-jet
pixi install
```

If you already have a copy of the repository, open that folder and run `pixi install` from there.

### Input data and references

- **Paired-end RNA-seq FASTQs** for each sample, with separate R1 and R2 paths.
- **HLA class I alleles** for each sample before running pVACbind. Separate alleles with commas, for example `HLA-A*11:01,HLA-B*13:01,HLA-C*07:02`.
- **A GRCh38 gene annotation GTF** and a matching **RepeatMasker table**. The RepeatMasker table describes the genome and can be shared by samples using the same assembly.
- **A STAR index** compatible with your STAR version, or a matching genome FASTA and GTF to build one. When reusing an index, point to the directory that directly contains `Genome`, `SA`, `SAindex` and `genomeParameters.txt`.

### Software and containers

| Workflow | What to prepare |
|---|---|
| Nextflow | Pixi, Apptainer, and the `preproc.sif`, `jet-r-4.5.3.sif` and `pVACtools.sif` images. The R image needs the JET analysis packages and the hg38 genome package. |
| Bash | Pixi for STAR and samtools, a host R installation with the JET packages available to `Rscript`, Apptainer, and the `pVACtools.sif` image. |

The pVACtools image needs the local IEDB installation at `/opt/iedb` inside the container. NetMHCStabPan can still require internet access for stability predictions.

The commands below use the development server's example paths. Replace them with paths on your machine. For Nextflow, also check the image paths in [`conf/base.config`](NF_JET_identificatio_pipeline/conf/base.config) and the work directory and Apptainer binds in [`nextflow.config`](NF_JET_identificatio_pipeline/nextflow.config).

## 5. Run with Nextflow

Prepare a tab-separated manifest with the header `sampleName`, `rnaFastq1`, `rnaFastq2`, `hlaAlleles`. Add one row per sample. For a multi-sample run, put all samples in the same manifest and pass that file once with `--manifestPath`. The GTF, STAR index and RepeatMasker table are supplied once for the run.

From the repository root, run:

```bash
data=/home/faramir/data/neondisco-jet
manifest="$PWD/infoDir/379T_nf_manifest.tsv"
out="$data/outputData/379T_nextflow_steps1to3"

pixi run nextflow run NF_JET_identificatio_pipeline/main.nf \
  -profile apptainer \
  --step all \
  --manifestPath "$manifest" \
  --starIndex "$data/starIndexes_2.7.11b/star_index" \
  --gtf "$PWD/inputData/Homo_sapiens.GRCh38.115.gtf" \
  --repeatsFile "$PWD/metadataDir/379T.repeatmasker_reformat.txt" \
  --genome hg38 \
  --outputDir "$out"
```

To run a cohort, change `manifest` to your multi-sample TSV and `out` to a new output folder. To build an index within the workflow, supply `--genomeFa` in place of `--starIndex`.

See the [Nextflow guide](NF_JET_identificatio_pipeline/README.md) for a two-sample manifest, individual step commands, resource settings and DAG creation.

## 6. Run with Bash

Prepare `infoDir/samples.tsv` with **five tab-separated columns and no header**. Each row contains the sample name, R1 path, R2 path, HLA alleles and RepeatMasker path, in that order. Use absolute paths in the sample sheet.

From the repository root, run:

```bash
repo="$PWD"
data=/home/faramir/data/neondisco-jet
pixi_bin="$repo/.pixi/envs/default/bin"

bash JET_identification_pipeline/run_TEAL_pipeline.sh \
  --samples "$repo/infoDir/samples.tsv" \
  --genome hg38 \
  --genome-fasta "$repo/inputData/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --gtf "$repo/inputData/Homo_sapiens.GRCh38.115.gtf" \
  --star-index "$data/starIndexesDir" \
  --star-bin "$pixi_bin" \
  --samtools-bin "$pixi_bin" \
  --rscript-bin /opt/R/4.5.3/bin/Rscript \
  --r-script "$repo/JET_identification_pipeline/JET_analysis_filtered.R" \
  --pvactools-sif /home/faramir/images/pVACtools.sif \
  --outputs-dir "$data/outputData" \
  --tmp-dir "$data/tmpData" \
  --threads 16
```

See the [Bash guide](JET_identification_pipeline/README.md) for R package setup, sample sheet examples, skip options and how to interpret the final combined table.

## 7. Outputs and Run Reports

| Workflow | Where to look |
|---|---|
| **Nextflow** | Under `--outputDir`, check `STAR-alignments/`, `JET-identification/` and `pVACbind/` for results by sample. The `reports/` folder contains the execution report, timeline and task trace. |
| **Bash** | Under `--outputs-dir`, check the dated sample folders for STAR outputs, JET tables, pVACbind results and `<sample>_final_combined_data_TEAL.tsv`. |

Keep the peptide ID tables alongside the pVACbind results. The ID tables let you trace a prediction back to its JET records.

If a Nextflow task fails, check the Nextflow log and `reports/trace.txt` for the process name and sample, then inspect `.command.err` in the failed task's work directory. Keep both the cache and work directory if you plan to rerun with `-resume`.

## 8. Citation and Credits

If you use this pipeline, please cite:

> Burbage M, Rocañín-Arjó A, et al. *Epigenetically-controlled tumor antigens derived from splice junctions between exons and transposable elements.* **Science Immunology**, 2023.

> Hundal J, et al. *pVACtools: a computational toolkit to identify and visualize cancer neoantigens.* **Cancer Immunology Research**, 2020.

**Original JET pipeline:** Alexandre Houy and Christel Goudot (Institut Curie), with modifications by Ares Rocañín-Arjó.

**GRCh38 adaptation, pVACbind integration and multi-sample workflow:** Muhammad Faramir (Cancer Research Malaysia).
