// enable dsl2
nextflow.enable.dsl=2

// import modules

workflow {

    // Set channel for the fastq directories
    barcodeDirs = file("${params.inputDir}/barcode*", type: 'dir', maxdepth: 1 )

    ch_sample = Channel
            .fromPath(barcodeDirs)
            .filter(~/.*barcode[0-9]{1,4}$/)
            .map{ dir ->
                def count = 0
                for (x in dir.listFiles()) {
                    if (x.isFile() && x.toString().contains('.fastq')) {
                        count += x.countFastq()
                    }
                }
                return[dir.baseName, dir, count]
            }




    main:

        ch_sample.view()

}