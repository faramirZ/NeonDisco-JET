process PREDICT_PVACBIND {
    tag "${sampleName}:size${peptideSize}"
    label 'predictPvacbind'

    container "${params.container__pvactools}"

    publishDir {
        "${params.outputDir}/pVACbind/${sampleName}"
    }, mode: 'copy', pattern: 'pvacbind_size*'

    input:
    tuple val(sampleName),
          val(peptideSize),
          path(fasta),
          path(ids),
          val(hlaAlleles)

    output:
    tuple val(sampleName),
          val(peptideSize),
          path("pvacbind_size${peptideSize}"),
          path(ids),
          emit: predictions

    script:
    if (!(peptideSize in [8, 9, 10, 11])) {
        error "Unsupported peptide size: ${peptideSize}"
    }

    """
    set -euo pipefail

    pvacbind-JET.sh \\
        "${sampleName}" \\
        "${fasta}" \\
        "${hlaAlleles}" \\
        ${peptideSize} \\
        ${task.cpus} \\
        "pvacbind_size${peptideSize}"
    """
}
