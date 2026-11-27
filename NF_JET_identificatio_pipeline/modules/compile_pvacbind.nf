process COMPILE_PVACBIND {
    tag "${sampleName}"
    cache false

    publishDir {
        "${params.outputDir}/finalresult/${sampleName}"
    }, mode: 'copy'

    input:
    tuple val(sampleName), val(resultDirs)

    output:
    tuple val(sampleName), path("finalresult${sampleName}.tsv"), emit: results

    exec:
    def reports = []
    def allColumns = []

    // Read each header first so the final table includes every column.
    resultDirs.each { dir ->
        def size = dir.name - 'pvacbind_size'
        def report = dir.resolve("MHC_Class_I/${sampleName}.MHC_I.filtered.tsv")
        if (!report.exists()) return

        def header = report.toFile().withReader('UTF-8') { reader ->
            reader.readLine()
        }
        if (!header) error "Missing header in ${report}"

        def columns = header.split('\t', -1).toList()
        if (columns.contains('') || columns.unique(false).size() != columns.size()) {
            error "Blank or duplicate column name in ${report}"
        }
        if (columns.contains('Size')) {
            error "Report already contains a column named Size: ${report}"
        }

        columns.each { column ->
            if (!allColumns.contains(column)) allColumns.add(column)
        }
        reports.add([size: size, file: report, columns: columns])
    }

    if (reports.isEmpty()) error "No filtered reports found for ${sampleName}"

    def outputFile = task.workDir.resolve("finalresult${sampleName}.tsv")
    outputFile.toFile().withWriter('UTF-8') { writer ->
        writer.write("Size\t${allColumns.join('\t')}\n")

        reports.each { report ->
            // Find each output column's position in this particular report.
            def positions = allColumns.collect { column ->
                report.columns.indexOf(column)
            }

            report.file.toFile().withReader('UTF-8') { reader ->
                reader.readLine()  // skip this report's header
                def lineNumber = 1

                reader.eachLine { line ->
                    lineNumber = lineNumber + 1
                    if (line.isEmpty()) return

                    def fields = line.split('\t', -1)
                    if (fields.length != report.columns.size()) {
                        error "Wrong number of fields in ${report.file} at line ${lineNumber}"
                    }

                    def values = positions.collect { position ->
                        position < 0 ? 'NA' : fields[position]
                    }
                    writer.write("${report.size}\t${values.join('\t')}\n")
                }
            }
        }
    }
}
