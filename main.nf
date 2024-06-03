// enable dsl2
nextflow.enable.dsl=2

// import modules

workflow {

    // Set channel for the fastq directories
    nanoporeBarcodeDirs = file("${params.inputDir}/barcode*", type: 'dir', maxdepth: 1 )

    Channel
            .fromPath(nanoporeBarcodeDirs)
            .filter(~/.*barcode[0-9]{1,4}$/)
            .filter{ d ->
                def count = 0
                for (x in d.listFiles()) {
                    if (x.isFile()) {
                        count += x.countFastq()
                    }
                }
                count > params.OntMinReadsPerBarcode
            }
            .set{ch_sample}
            .ifEmpty{ error "Cannot find any barcode directories matching: ${params.inputDir}"
            }



    main:

        ch_sample.view()

}