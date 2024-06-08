// enable dsl2
nextflow.enable.dsl=2

// import modules
include {concatenate} from './modules/concatenate.nf'
include {nanoq} from './modules/nanoq.nf'
include {fastqcRaw} from './modules/fastqc.nf'
include {fastqcTrimmed} from './modules/fastqc.nf'
include {minimap} from './modules/minimap.nf'
include {ivar} from './modules/ivar.nf'
include {sortIndex} from './modules/sortIndex.nf'
include {medaka} from './modules/medaka.nf'


workflow {

    // Set channel for the fastq directories
    barcodeDirs = file("${params.inputDir}/barcode*", type: 'dir', maxdepth: 1 )
    fastqDir = file("${params.inputDir}/*.fastq*" , type: 'file', maxdepth: 1)

    main:
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
                        return[dir.baseName, dir]
                    }
            
           
            concatenate(ch_sample)
            fastqcRaw(concatenate.out.concatFastq)
            nanoq(concatenate.out.concatFastq)


        } else if (fastqDir) {
            ch_sample = Channel
                    .fromPath(fastqDir)
                    .filter(file -> file.name =~ /.*\.fastq(\.gz)?$/)
                    .map{file ->
                        def baseName = file.name.replaceAll(/\.fastq(\.gz)?$/, '')
                        def reads = file.countFastq()
                        return [baseName, file]
                    }
            fastqcRaw(ch_sample)
            nanoq(ch_sample)

        } else {
            log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files or the concatenated fastq files e.g. --inputDir ./raw/fastq_pass/ or --inputDir ./fastqConcatenated/"
            System.exit(1)
        }

        fastqcTrimmed(nanoq.out.trimmedFastq)
        minimap(nanoq.out.trimmedFastq, params.reference)
        ivar(minimap.out.bam)
        sortIndex(ivar.out.trimmedBam)
        medaka(sortIndex.out.bamBai)




}