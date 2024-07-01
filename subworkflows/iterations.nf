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

include {minimap2Iterate} from '../modules/ont/minimap2.nf'




workflow stage1 {
    take:
        ch_trimmedFastq

    main:
        minimap2(ch_trimmedFastq, params.reference)
        sam2bam(minimap2.out.sam)
        sortIndexMinimap(sam2bam.out.bam)
        trimPrimer(sortIndexMinimap.out.bamBai)
        sortIndexIvar(trimPrimer.out.trimmedBam)
        medaka(sortIndexIvar.out.bamBai)

    emit:
        consensus           =           medaka.out.consensus
}


workflow stage2 {
    take:
        ch_trimmedFastq
        ch_consensus

    main:

        minimap2Iterate(ch_trimmedFastq.join(ch_consensus))
        sam2bam(minimap2Iterate.out.sam)
        sortIndexMinimap(sam2bam.out.bam)
        trimPrimer(sortIndexMinimap.out.bamBai)
        sortIndexIvar(trimPrimer.out.trimmedBam)
        medaka(sortIndexIvar.out.bamBai)

    emit:
        consensus           =           medaka.out.consensus
}


workflow stage3 {
    take:
        ch_trimmedFastq
        ch_consensus

    main:

        minimap2Iterate(ch_trimmedFastq.join(ch_consensus))
        sam2bam(minimap2Iterate.out.sam)
        sortIndexMinimap(sam2bam.out.bam)
        trimPrimer(sortIndexMinimap.out.bamBai)
        sortIndexIvar(trimPrimer.out.trimmedBam)
        medaka(sortIndexIvar.out.bamBai)

    emit:
        consensus           =           medaka.out.consensus
}
