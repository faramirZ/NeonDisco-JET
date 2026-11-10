process BUILD_STAR_INDEX {
    tag "${indexName}"
    label 'buildStarIndex'

    container "${params.container__preproc}"
    publishDir "${params.starIndexesDir}", mode: 'copy'

    input:
    path genomeFa
    path gtf
    val indexName
    val sjdbOverhang

    output:
    path "${indexName}", emit: index

    script:
    """
    set -euo pipefail

    FA="${genomeFa}"
    GTF="${gtf}"

    if [[ "\$FA" == *.gz ]]; then
        gunzip -c "\$FA" > genome.fa
        FA=genome.fa
    fi

    if [[ "\$GTF" == *.gz ]]; then
        gunzip -c "\$GTF" > annotation.gtf
        GTF=annotation.gtf
    fi

    mkdir -p "${indexName}"

    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir "${indexName}" \\
        --genomeFastaFiles "\$FA" \\
        --sjdbGTFfile "\$GTF" \\
        --sjdbOverhang ${sjdbOverhang}
    """
}
