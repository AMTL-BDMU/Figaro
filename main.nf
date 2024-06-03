// enable dsl2
nextflow.enable.dsl=2

// import modules

workflow {

    // Set channel for the fastq directories
    nanoporeBarcodeDirs = file("${params.inputDir}/barcode*", type: 'dir', maxdepth: 1 )

    ch_sample = Channel
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
            .ifEmpty{exit 1, "Cannot find any barcode directories."
            }



    main:

        ch_sample.view()

}