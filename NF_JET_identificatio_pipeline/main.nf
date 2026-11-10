nextflow.enable.dsl = 2

include { BUILD_STAR_INDEX } from './modules/build_star_index'
include { ALIGN_READS_JET } from './modules/align_reads_jet'
include { IDENTIFY_JETS } from './modules/identify_jets'
include { PREDICT_PVACBIND } from './modules/predict_pvacbind'

params.step = 'all'

workflow {
    def step = params.step.toString()

    if (!(step in ['all', '1', '2', '3'])) {
        error '--step must be all, 1, 2, or 3'
    }

    if (!params.manifestPath) {
        error 'Provide --manifestPath'
    }

    // Validate the manifest before launching sample tasks.
    samples = channel
        .fromPath(params.manifestPath, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .toList()
        .flatMap { rows ->
            if (!rows) error 'The manifest contains no samples'

            def seen = new HashSet()

            rows.each { row ->
                if (!(row.sampleName ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                    error "Invalid sampleName: ${row.sampleName}"
                }

                if (!seen.add(row.sampleName)) {
                    error "Duplicate sample: ${row.sampleName}"
                }

                if (step in ['all', '3'] && !row.hlaAlleles?.trim()) {
                    error "Missing hlaAlleles for ${row.sampleName}"
                }
            }

            rows
        }

    if (step in ['all', '1', '2']) {
        if (!params.gtf) error 'Provide --gtf'
        gtf = file(params.gtf, checkIfExists: true)
    }

    // STEP 1: shared index and independent sample alignments.
    if (step in ['all', '1']) {
        reads = samples.map { row ->
            if (!row.rnaFastq1 || !row.rnaFastq2) {
                error "Missing FASTQ paths for ${row.sampleName}"
            }

            tuple(row.sampleName, [
                file(row.rnaFastq1, checkIfExists: true),
                file(row.rnaFastq2, checkIfExists: true)
            ])
        }

        if (params.starIndex) {
            def index = file(params.starIndex, checkIfExists: true)

            ['Genome', 'SA', 'SAindex', 'genomeParameters.txt'].each { name ->
                if (!index.resolve(name).exists()) {
                    error "Index missing ${name}: ${index}"
                }
            }

            index_ch = channel.value(index)
        } else {
            if (!params.genomeFa) {
                error 'Provide --starIndex or --genomeFa'
            }

            BUILD_STAR_INDEX(
                file(params.genomeFa, checkIfExists: true),
                gtf,
                params.starIndexName,
                params.sjdbOverhang
            )

            index_ch = BUILD_STAR_INDEX.out.index
        }

        ALIGN_READS_JET(reads, index_ch, gtf, params.sjdbOverhang)
        star_outputs = ALIGN_READS_JET.out.jet_inputs
    }

    // Standalone Step 2: load existing STAR outputs.
    if (step == '2') {
        star_outputs = samples.map { row ->
            if (!row.starDir) {
                error "Missing starDir for ${row.sampleName}"
            }

            def prefix = "${row.starDir}/${row.sampleName}"

            tuple(
                row.sampleName,
                file("${prefix}_Chimeric.out.junction", checkIfExists: true),
                file("${prefix}_SJ.out.tab", checkIfExists: true),
                file("${prefix}_Log.final.out", checkIfExists: true)
            )
        }
    }

    // STEP 2: process each sample as its STAR outputs arrive.
    if (step in ['all', '2']) {
        if (!params.repeatsFile || !params.jetRScript || !params.genome) {
            error 'Provide --repeatsFile, --jetRScript, and --genome'
        }

        IDENTIFY_JETS(
            star_outputs,
            file(params.repeatsFile, checkIfExists: true),
            gtf,
            file(params.jetRScript, checkIfExists: true),
            params.genome
        )

        peptide_outputs = IDENTIFY_JETS.out.peptides
    }

    // Standalone Step 3: load existing FASTA and ID files.
    if (step == '3') {
        peptide_outputs = samples.map { row ->
            if (!row.jetDir) {
                error "Missing jetDir for ${row.sampleName}"
            }

            def prefix =
                "${row.jetDir}/${row.sampleName}_Fusions_chim2.junc2.size"

            def fastas = (8..11).collect { size ->
                file("${prefix}${size}.fasta", checkIfExists: true)
            }

            def ids = (8..11).collect { size ->
                file("${prefix}${size}.ids.txt", checkIfExists: true)
            }

            tuple(row.sampleName, fastas, ids)
        }
    }

    // STEP 3: match HLA by sample, then match FASTA and IDs by size.
    if (step in ['all', '3']) {
        hla = samples.map { row ->
            tuple(row.sampleName, row.hlaAlleles)
        }

        prediction_inputs = peptide_outputs
            .join(hla, failOnDuplicate: true, failOnMismatch: true)
            .flatMap { sampleName, fastas, ids, alleles ->
                (8..11).collect { size ->
                    def stem =
                        "${sampleName}_Fusions_chim2.junc2.size${size}"

                    def fasta = fastas.find {
                        it.name == "${stem}.fasta"
                    }

                    def idFile = ids.find {
                        it.name == "${stem}.ids.txt"
                    }

                    if (!fasta || !idFile) {
                        error "Missing peptide pair: ${sampleName}, size ${size}"
                    }

                    tuple(sampleName, size, fasta, idFile, alleles)
                }
            }

        PREDICT_PVACBIND(prediction_inputs)
    }
}
