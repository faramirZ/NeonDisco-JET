process QUANT_COUNTS_KALLISTO {
    tag "${sampleName}"
	label 'quantCounts'
    container "${params.container__kallisto}"
    publishDir "${params.outputDir}/kallisto_gene_counts", mode: 'copy'
    
    cpus 8
    
    input:
    tuple val(sampleName), path(read1), path(read2)
    path(index)
    
    output:
    path("${sampleName}/${sampleName}_abundance.tsv"), emit: abundance
    path("${sampleName}/${sampleName}_abundance.h5"), emit: h5
    path("${sampleName}/${sampleName}_run_info.json"), emit: run_info
    
    script:
    """
    if kallisto quant \
        -i ${index} \
        -o ${sampleName} \
        -b ${params.bootstrap_samples} \
        -t ${task.cpus} \
        ${read1} ${read2};
    then
        echo "Kallisto quantification completed successfully for sample ${sampleName}"
        # rename output files to include sample name
        mv ${sampleName}/abundance.tsv ${sampleName}/${sampleName}_abundance.tsv && \
        mv ${sampleName}/abundance.h5 ${sampleName}/${sampleName}_abundance.h5 && \
        mv ${sampleName}/run_info.json ${sampleName}/${sampleName}_run_info.json
    else
        echo "Kallisto quantification failed for sample ${sampleName}" >&2
        exit 1
    fi

    """
}
