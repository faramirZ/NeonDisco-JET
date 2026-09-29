process QC_READS_FASTQC {
    tag "$sampleId (${stage})"

    label 'qcReadsFastQC'
    container "${params.container__qc}"

    publishDir { "${params.outputDir}/fastqc/${stage}" }, mode: 'copy', pattern: "*_fastqc.*"

    input:
    tuple val(sampleId), path(read1), path(read2)
    val stage // "raw" or "trimmed"

    output:
    tuple val(sampleId), path("*_fastqc.zip"), path("*_fastqc.html"), emit: fastqc_files

    script:
    """
    fastqc ${read1} ${read2}
    """
}
