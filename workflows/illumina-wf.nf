// enable dsl2
nextflow.enable.dsl=2


// import modules
include {fastP} from '../modules/illumina/fastP.nf'
include {fastqcRawPE} from '../modules/misc/fastqc.nf'
include {fastqcTrimmedPE} from '../modules/misc/fastqc.nf'
include {hydra} from '../modules/illumina/hydra.nf'
include {sierra} from '../modules/misc/sierra.nf'
include {pdfReport} from '../modules/misc/pdfReport.nf'

workflow illuminaShotgun {
    Channel
        .fromFilePairs("${params.inDir}/*{,.trimmed}_{R1,R2,1,2}{,_001}.{fastq,fq}{,.gz}", flat:true)
        .ifEmpty{error "Cannot find any reads matching: ${params.inDir}"}
        .set{ch_sample}
        

    main:
        fastqcRawPE(ch_sample)
        fastP(ch_sample)
        fastqcTrimmedPE(fastP.out.trimmed)
        hydra(fastP.out.trimmed)
        sierra(hydra.out.consensus)
        pdfReport(sierra.out.json, params.markdownFile)

}