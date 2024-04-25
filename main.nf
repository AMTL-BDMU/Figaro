// enable dsl2
nextflow.enable.dsl=2

// import modules
include { cat } from './modules/01-nanoq.nf'
include { nanoq } from './modules/01-nanoq.nf'

include { minimap2 } from './modules/02-reference_map.nf'
include { samtools } from './modules/02-reference_map.nf'
include { nanopolish } from './modules/02-reference_map.nf'

include { minimap2_2 } from './modules/03-consensus_map.nf'
include { samtools_2 } from './modules/03-consensus_map.nf'
include { nanopolish_2 } from './modules/03-consensus_map.nf'

include { minimap2_3 } from './modules/04-consensus_map.nf'
include { samtools_3 } from './modules/04-consensus_map.nf'
include { nanopolish_3 } from './modules/04-consensus_map.nf'

include { bgzip } from './modules/04-consensus_map.nf'

include { consens_bcf } from './modules/05-consensus_mixedSite.nf'
include { consens_quasi } from './modules/05-consensus_mixedSite.nf'
include { consens_varscan } from './modules/05-consensus_mixedSite.nf'
include { consens_lofreq } from './modules/05-consensus_mixedSite.nf'

include { sierra } from './modules/06-sierra.nf' 
include { reportDrugResistance } from './modules/07-reportDrugResistance.nf'
include { drugResistanceScoreCombine } from './modules/08-drugResistanceScoreCombine.nf'



workflow {

         ch_sample = Channel
                .fromPath("${params.reads}/**/fastq_pass/barcode*/", type: 'dir')
                .ifEmpty { error "Cannot find any barcode directories matching: ${params.reads}" }

        main:
                ch_sample.map { barcodePath -> tuple(barcodePath) }

                cat( ch_sample )
                nanoq( cat.out.cat_fastq )

                minimap2( nanoq.out.nanoq_fastq, params.reference)
                samtools( minimap2.out.ref_aligned_sam)
                nanopolish( nanoq.out.nanoq_fastq.join(samtools.out.ref_aligned_bam), params.reference)

                minimap2_2( nanopolish.out.first_consensus.join(nanoq.out.nanoq_fastq) )
                samtools_2( minimap2_2.out.con1_aligned_sam )
                nanopolish_2( (((nanoq.out.nanoq_fastq).join(nanopolish.out.first_consensus)).join(samtools_2.out.con1_aligned_bam)).view(), params.reference )

                minimap2_3( nanopolish_2.out.second_consensus.join(nanoq.out.nanoq_fastq) )
                samtools_3( minimap2_3.out.con2_aligned_sam )
                nanopolish_3( ((((nanoq.out.nanoq_fastq).join(nanopolish.out.first_consensus)).join(samtools_3.out.con2_aligned_bam))).join(nanopolish_2.out.second_consensus), params.reference )

                bgzip ( nanopolish_3.out.final_consensus )
                consens_bcf( (nanopolish_3.out.final_consensus).join(bgzip.out.bgzip_vcf) )
                consens_quasi( (nanopolish_3.out.final_consensus).join(bgzip.out.bgzip_vcf) )
                consens_varscan( (nanopolish_3.out.final_consensus).join(bgzip.out.bgzip_vcf) )
                consens_lofreq( (nanopolish_3.out.final_consensus).join(bgzip.out.bgzip_vcf) )

                sierra( nanopolish_3.out.final_consensus )

                drugResistanceScoreCombine(sierra.out.csv.collect())

                reportDrugResistance(sierra.out.json, params.reportPDF)

}
