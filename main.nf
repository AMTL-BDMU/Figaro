// enable dsl2
nextflow.enable.dsl=2

// import modules
include { cat } from './modules/nanoq.nf'
include { nanoq } from './modules/nanoq.nf'

include { downloadRef } from './modules/reference_map.nf'
include { minimap2 } from './modules/reference_map.nf'
include { samtools } from './modules/reference_map.nf'
include { medaka } from './modules/reference_map.nf'

include { minimap2_2 } from './modules/consensus_map.nf'
include { samtools_2 } from './modules/consensus_map.nf'
include { medaka_2 } from './modules/consensus_map.nf'

include { sierra } from './modules/sierra.nf'
 
include {reportDrugResistance} from './modules/reportDrugResistance.nf'
include {drugResistanceScoreCombine} from './modules/drugResistanceScoreCombine.nf'

workflow {

        Channel

                .fromPath("${params.reads}/**/fastq_pass/barcode*/", type: 'dir')
                .ifEmpty { error "Cannot find any barcode directories matching: ${params.reads}" }
                .set { ch_sample }


        main:
                ch_sample.map { barcodePath -> tuple(barcodePath) }

                cat( ch_sample )
                nanoq( cat.out.cat_fastq )

                downloadRef()

                minimap2( nanoq.out.nanoq_fastq, downloadRef.out.reference )
                samtools( minimap2.out.ref_aligned_sam, downloadRef.out.reference)
                medaka( samtools.out.first_fasta, downloadRef.out.reference )

                minimap2_2( medaka.out.first_fasta_polished, nanoq.out.nanoq_fastq )
                samtools_2( minimap2_2.out.con1_aligned_sam, medaka.out.first_fasta_polished )
                medaka_2( samtools_2.out.second_fasta, medaka.out.first_fasta_polished )

                sierra( medaka_2.out.second_fasta_polished )
                drugResistanceScoreCombine(sierra.out.csv.collect())
                reportDrugResistance(sierra.out.json, params.reportPDF)

}
