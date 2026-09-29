process BUILD_STAR_INDEX {
    tag "${indexName}"
    label 'buildStarIndex'
    container "${params.container__preproc}"
    publishDir "${params.starIndexesDir}", mode: 'copy'

    input:
    path genomeFa
    val  gtf          // plain path string, or '' for no annotation
    val  indexName
    val  sjdbOverhang

    output:
    path indexName, emit: index

    script:
    """
    set -euo pipefail

    # STAR needs plain (uncompressed) FASTA/GTF
    FA="${genomeFa}"
    if [[ "\$FA" == *.gz ]]; then
        gunzip -c "\$FA" > genome.fa
        FA=genome.fa
    fi

    GTF_ARGS=""
    if [[ -n "${gtf}" ]]; then
        GTF="${gtf}"
        if [[ "\$GTF" == *.gz ]]; then
            gunzip -c "\$GTF" > annotation.gtf
            GTF=annotation.gtf
        fi
        GTF_ARGS="--sjdbGTFfile \$GTF --sjdbOverhang ${sjdbOverhang}"
    fi

    mkdir -p ${indexName}

    STAR --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir ${indexName} \
        --genomeFastaFiles "\$FA" \
        \$GTF_ARGS
    """
}
