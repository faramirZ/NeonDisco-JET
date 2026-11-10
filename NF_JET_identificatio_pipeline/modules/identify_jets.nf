process IDENTIFY_JETS {
    tag "${sampleName}"
    label 'identifyJets'

    publishDir {
        "${params.outputDir}/JET-identification/${sampleName}"
    }, mode: 'copy'

    input:
    tuple val(sampleName),
          path(chimeric),
          path(junctions),
          path(starLog)

    path repeats
    path gtf
    path rScript
    val genome

    output:
    tuple val(sampleName),
          path("${sampleName}_Chimeric.out.annotatedJET.txt"),
          emit: annotated_jets

    tuple val(sampleName),
          path("${sampleName}_Fusions_chim2.junc2.size*.fasta",
               arity: '4'),
          path("${sampleName}_Fusions_chim2.junc2.size*.ids.txt",
               arity: '4'),
          emit: peptides

    tuple val(sampleName),
          path("${sampleName}_Fusions_chim2.junc2.size*.genomic.txt",
               arity: '4'),
          emit: genomic_sequences

    tuple val(sampleName),
          path("${sampleName}_Fusions_chim2.junc2.size*.RData",
               arity: '4'),
          emit: rdata

    script:
    """
    set -euo pipefail

    jet-identify.sh \\
        "${sampleName}" \\
        "${chimeric}" \\
        "${junctions}" \\
        "${starLog}" \\
        "${repeats}" \\
        "${gtf}" \\
        "${rScript}" \\
        "${genome}"
    """
}
