#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Simplified pipeline: raw QC -> fastp trimming -> trimmed QC -> STAR alignment
 * Containers:
 *   - preproc.sif  -> FastQC + fastp                 (params.container__preproc)
 *   - star.sif     -> STAR index build + alignment   (params.container__star)
 *
 * STAR index: searched for automatically; built only if none is found.
 */

include { QC_READS_FASTQC  as QC_RAW     } from './modules/qc_reads_fastqc.nf'
include { QC_READS_FASTQC  as QC_TRIMMED } from './modules/qc_reads_fastqc.nf'
include { TRIM_READS_FASTP               } from './modules/trim_reads_fastp.nf'
include { BUILD_STAR_INDEX               } from './modules/build_star_index.nf'
include { ALIGN_READS_STAR               } from './modules/align_reads_star.nf'

// A directory counts as a STAR index if the core index files are present
def isValidStarIndex(File d) {
    d.isDirectory() && ['genomeParameters.txt', 'Genome', 'SA', 'SAindex'].every { new File(d, it).exists() }
}

// Search order: --starIndex, then starIndexesDir itself, then each subfolder of starIndexesDir
def findStarIndex(String indexesDir, String explicitIndex) {
    def candidates = []
    if (explicitIndex) candidates << new File(explicitIndex)
    def root = new File(indexesDir)
    candidates << root
    if (root.isDirectory()) {
        (root.listFiles() ?: []).findAll { it.isDirectory() }.sort { it.name }.each { candidates << it }
    }
    return candidates.find { isValidStarIndex(it) }
}

workflow {

    if (!params.manifestPath) {
        error "--manifestPath must be provided"
    }

    // ---- FASTQ channel from a TSV manifest ----
    // Expected columns: sampleName  rnaFastq1  rnaFastq2
    fastq_ch = channel
        .fromPath(params.manifestPath)
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.rnaFastq1 && row.rnaFastq2 }
        .map { row ->
            def read1 = file(row.rnaFastq1)
            def read2 = file(row.rnaFastq2)

            // allows e.g. 379T_R1_trim.fastq.gz as well as sample_R1.fastq.gz
            def r1ok = read1.name =~ /[_\.][Rr]?1(_[A-Za-z]+)?\.(fastq|fq)(\.gz)?$/
            def r2ok = read2.name =~ /[_\.][Rr]?2(_[A-Za-z]+)?\.(fastq|fq)(\.gz)?$/
            if (!r1ok || !r2ok) {
                error "Invalid read pairing format for sample ${row.sampleName}: ${read1.name}, ${read2.name}"
            }
            tuple(row.sampleName, read1, read2)
        }

    fastq_ch.view { "FASTQ input: $it" }

    // ---- STAR index: find it, or build it ----
    def foundIndex = findStarIndex(params.starIndexesDir, params.starIndex)

    if (foundIndex) {
        log.info "Using existing STAR index: ${foundIndex}"
        star_index_ch = channel.value(file(foundIndex.toString()))
    }
    else {
        if (!params.genomeFa) {
            error "No STAR index found in ${params.starIndexesDir}. Provide --genomeFa (and ideally --gtf) so one can be built."
        }
        log.info "No STAR index found in ${params.starIndexesDir}; building '${params.starIndexName}' from ${params.genomeFa}"
        BUILD_STAR_INDEX(
            file(params.genomeFa),
            params.gtf ? params.gtf.toString() : '',
            params.starIndexName,
            params.sjdbOverhang
        )
        star_index_ch = BUILD_STAR_INDEX.out.index
    }

    // ---- QC on raw reads ----
    QC_RAW(fastq_ch, "raw")

    // ---- Trimming ----
    TRIM_READS_FASTP(fastq_ch)

    // ---- QC on trimmed reads ----
    QC_TRIMMED(TRIM_READS_FASTP.out.trimmed_reads, "trimmed")

    // ---- STAR alignment ----
    // Reshape tuple(sampleName, r1, r2) -> tuple(sampleName, [r1, r2])
    star_input_ch = TRIM_READS_FASTP.out.trimmed_reads
        .map { sampleName, r1, r2 -> tuple(sampleName, [r1, r2]) }

    ALIGN_READS_STAR(star_input_ch, star_index_ch)

    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Duration: $workflow.duration"
    }
}
