process ALIGN_READS_JET {
    tag "${sampleName}"
    label 'alignReadsStar'

    container "${params.container__preproc}"
    publishDir {
        "${params.outputDir}/STAR-alignments/${sampleName}/processed_calls"
    }, mode: 'copy'

    input:
    tuple val(sampleName), path(reads, stageAs: 'reads??/*')
    path starIndex
    path gtf
    val sjdbOverhang

    output:
    tuple val(sampleName),
          path("${sampleName}_Chimeric.out.junction"),
          path("${sampleName}_SJ.out.tab"),
          path("${sampleName}_Log.final.out"),
          emit: jet_inputs

    tuple val(sampleName),
          path("${sampleName}_Aligned.sortedByCoord.out.bam"),
          path("${sampleName}_Aligned.sortedByCoord.out.bam.bai"),
          emit: final_bam

    tuple val(sampleName),
          path("${sampleName}_ReadsPerGene.out.tab"),
          emit: gene_counts

    script:
    if (reads.size() != 2) {
        error "Expected two FASTQ files for ${sampleName}"
    }

    """
    set -euo pipefail

    star-align-JET.sh \\
        "${reads[0]}" \\
        "${reads[1]}" \\
        "${sampleName}" \\
        ${task.cpus} \\
        "${starIndex}" \\
        "${gtf}" \\
        ${sjdbOverhang}

    samtools index "${sampleName}_Aligned.sortedByCoord.out.bam"
    """
}
