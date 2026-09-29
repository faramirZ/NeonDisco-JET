// Defines the FastQC process
process QC_READS_FASTQC {
    tag "$sample_id"
    label 'qcReadsFastQC'
    container "${params.container__kallisto}"
    publishDir "${params.outputDir}/fastqc/${sample_id}", mode: 'copy', pattern: "*fastqc.*"

    input:
    // Receives the existing tuple structure
    tuple val(sample_id), path(read1), path(read2)

    output:
    // Output the HTML/ZIP files
    tuple val(sample_id), path("*fastqc.zip"), path("*fastqc.html")

    script:
    """
    # fastqc accepts R1 and R2 files separated by a space
    fastqc ${read1} ${read2}
    """
}