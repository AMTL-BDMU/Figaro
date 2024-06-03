// enable dsl2
nextflow.enable.dsl=2

// import modules

workflow {

    // Set channel for the fastq directories
    barcodeDirs = file("${params.inputDir}/barcode*", type: 'dir', maxdepth: 1 )
    fastqDir = file("${params.inputDir}/*.fastq*" , type: 'file', maxdepth: 1)




    if (barcodeDirs) {
        ch_sample = Channel
                .fromPath(barcodeDirs)
                .filter(~/.*barcode[0-9]{1,4}$/)
                .map{dir ->
                    def reads = 0
                    for (x in dir.listFiles()) {
                        if (x.isFile() && x.toString().contains('.fastq')) {
                            reads += x.countFastq()
                        }
                    }
                    return[dir.baseName, dir, reads]
                }
    } else if (fastqDir) {
        ch_sample = Channel
                .fromPath("${params.inputDir}", type: 'dir', maxDepth: 1)
                .map{it -> [ 'SAMPLE_1', 'single_barcode', it, 10000000 ]}
    } else {
        log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files  or the concatenated fastq files \
        e.g. '--inputDir ./raw/fastq_pass/ or ./fastqConcatenated/"
        System.exit(1)
    }

    main:

        ch_sample.view()

}