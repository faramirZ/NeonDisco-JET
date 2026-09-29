process TRIM_READS_FASTP {
    tag "${sampleName}"
    label 'trimReads'
    container "${params.container__kallisto}"
    publishDir "${params.outputDir}/trimmed_reads", mode: 'copy', pattern: "*.html"
    publishDir "${params.outputDir}/trimmed_reads", mode: 'copy', pattern: "*_trimmed.fastq.gz"
    
    cpus 4
    
    input:
    tuple val(sampleName), path(read1), path(read2)
    
    output:
    tuple val(sampleName), path("${sampleName}_R1_trimmed.fastq.gz"), path("${sampleName}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    path("${sampleName}_fastp.html"), emit: html_report
    
    script:
    """
    fastp \
        -i ${read1} \
        -I ${read2} \
        -o ${sampleName}_R1_trimmed.fastq.gz \
        -O ${sampleName}_R2_trimmed.fastq.gz \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
		--overrepresentation_analysis \
        --html ${sampleName}_fastp.html
    """
}
