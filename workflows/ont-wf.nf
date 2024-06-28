// enable dsl2
nextflow.enable.dsl=2

// import modules
include {concatenate} from '../modules/ont/concatenate.nf'
include {nanoq} from '../modules/ont/nanoq.nf'
include {fastqcRaw} from '../modules/ont/fastqc.nf'
include {fastqcTrimmed} from '../modules/ont/fastqc.nf'



include {minimap2} from '../modules/ont/minimap2.nf'
include {sam2bam} from '../modules/ont/samtools.nf'
include {sortIndexMinimap} from '../modules/ont/samtools.nf'
include {trimPrimer} from '../modules/ont/ivar.nf'
include {sortIndexIvar} from '../modules/ont/samtools.nf'
include {medaka} from '../modules/ont/medaka.nf'

// import subworkflow
include {stage1} from '../subworkflows/iterations.nf'
include {stage2} from '../subworkflows/iterations.nf'
include {stage3} from '../subworkflows/iterations.nf'


workflow ontAmplicon {
    // Set channel for the fastq directories
    barcodeDirs = file("${params.inDir}/barcode*", type: 'dir', maxdepth: 1 )
    fastqDir = file("${params.inDir}/*.fastq*" , type: 'file', maxdepth: 1)

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
            fastqcTrimmed(nanoq.out.trimmedFastq)


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
            fastqcTrimmed(nanoq.out.trimmedFastq)

        } else {
            log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files or the concatenated fastq files e.g. --inDir ./raw/fastq_pass/ or --inDir ./fastqConcatenated/"
            System.exit(1)
        }

        stage1(nanoq.out.trimmedFastq)
        stage2(nanoq.out.trimmedFastq, stage1.out.consensus)
        stage3(nanoq.out.trimmedFastq, stage2.out.consensus)
}