// enable dsl2
nextflow.enable.dsl=2

// import modules
include { cat } from './modules/01-preprocess.nf'
include { nanoq } from './modules/01-preprocess.nf'
//include { fastqc } from './modules/01-preprocess.nf'
//include { multiqc } from './modules/01-preprocess.nf'

//include { medaka } from './modules/02-consensus.nf'

//include { sierra } from './modules/03-profiling.nf'

workflow {

         ch_sample = Channel
                .fromPath("${params.in_dir}/**/fastq_pass/barcode*/", type: 'dir')
                .ifEmpty { error "Cannot find any barcode directories matching: ${params.reads}" }
                

        main:
//               ch_sample.map { samplePath -> tuple(samplePath) }

                cat( ch_sample )
                nanoq( cat.out.cat_fastq )

}