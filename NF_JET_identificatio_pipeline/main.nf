nextflow.enable.dsl = 2

include { BUILD_STAR_INDEX } from './modules/build_star_index'
include { ALIGN_READS_JET } from './modules/align_reads_jet'
include { IDENTIFY_JETS } from './modules/identify_jets'
include { PREDICT_PVACBIND } from './modules/predict_pvacbind'
include { COMPILE_PVACBIND } from './modules/compile_pvacbind'

params.step = 'all'

workflow {
    def step = params.step.toString()

    if (!(step in ['all', '1', '2', '3', '4'])) {
        error '--step must be all, 1, 2, 3, or 4'
    }

    // Steps 1–3 read sample information from the manifest.
    // Standalone Step 4 reads the sample name from the filtered report filenames.
    if (step != '4') {
        if (!params.manifestPath) {
            error 'Provide --manifestPath'
        }

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
    }

    if (step in ['all', '1', '2']) {
        if (!params.gtf) error 'Provide --gtf'
        gtf = file(params.gtf, checkIfExists: true)
    }

    // STEP 1: use a shared STAR index and align each sample.
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

    // STEP 2: identify JETs as each sample's STAR outputs arrive.
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

    // Standalone Step 3: load existing peptide FASTAs and ID files.
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

    // STEP 3: match HLA by sample and run pVACbind once per peptide size.
    if (step in ['all', '3']) {
        hla = samples.map { row ->
            tuple(row.sampleName, row.hlaAlleles)
        }

        prediction_inputs = peptide_outputs
            .join(hla, failOnDuplicate: true, failOnMismatch: true)
            .flatMap { sampleName, fastas, ids, alleles ->
                (8..11).collect { size ->
                    def stem = "${sampleName}_Fusions_chim2.junc2.size${size}"

                    def fasta = fastas.find { path ->
                        path.name == "${stem}.fasta"
                    }

                    def idFile = ids.find { path ->
                        path.name == "${stem}.ids.txt"
                    }

                    if (!fasta || !idFile) {
                        error "Missing peptide pair: ${sampleName}, size ${size}"
                    }

                    tuple(sampleName, size, fasta, idFile, alleles)
                }
            }

        PREDICT_PVACBIND(prediction_inputs)

        // Group the four Step 3 outputs by sample. Each sample can move
        // to Step 4 as soon as its own four pVACbind tasks finish.
        results_by_sample = PREDICT_PVACBIND.out.predictions
            .map { sample, size, dir, _ids ->
                def report = dir.resolve("MHC_Class_I/${sample}.MHC_I.all_epitopes.tsv")
                def hasRows = report.exists() && report.toFile().withReader { reader ->
                    reader.readLine()          // skip the header
                    reader.readLine() != null  // at least one prediction
                }

                tuple(sample, tuple(size, dir), hasRows)
            }
            .groupTuple(size: 4)
            .map { sample, sizeDirs, results ->
                if (!results.any { hasRows -> hasRows }) {
                    error "pVACbind produced no predictions for ${sample} across sizes 8 to 11"
                }

                def dirs = sizeDirs
                    .sort { pair -> pair[0] }
                    .collect { pair -> pair[1] }

                tuple(sample, dirs)
            }
    }

    // Standalone Step 4: load existing pVACbind results for one sample.
    if (step == '4') {
        if (!params.pvacbindDir) {
            error 'Provide --pvacbindDir with the pVACbind results for one sample'
        }

        def sampleDir = file(params.pvacbindDir, checkIfExists: true)
        if (!sampleDir.isDirectory()) {
            error "Expected a pVACbind results directory: ${sampleDir}"
        }

        def dirs = (8..11)
            .collect { size -> sampleDir.resolve("pvacbind_size${size}") }
            .findAll { dir -> dir.isDirectory() }

        def suffix = '.MHC_I.filtered.tsv'
        def sampleNames = dirs
            .collectMany { dir -> files("${dir}/MHC_Class_I/*${suffix}") }
            .collect { report -> report.name.substring(0, report.name.size() - suffix.size()) }
            .unique()

        if (sampleNames.size() != 1 || !(sampleNames[0] ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Expected filtered pVACbind reports for exactly one sample in ${sampleDir}"
        }

        results_by_sample = channel.of(tuple(sampleNames[0], dirs))
    }

    // STEP 4: compile one final filtered TSV per sample.
    if (step in ['all', '4']) {
        COMPILE_PVACBIND(results_by_sample)
    }
}
