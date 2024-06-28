process ivar {
        container 'staphb/ivar:1.4.2'

        tag "Trimming primers out of ${sample}"

        publishDir (
        path: "${params.outDir}/04_ivarTrim",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(bam)

        output:
        tuple val(sample), path("*.primerTrimmed.bam"), emit:trimmedBam

        script:
        """
        ivar trim \
            -b $params.primer \
            -p ${sample}.primerTrimmed \
            -i ${bam} \
            -q 1 \
            -s 4
        """
}


ivar trim -i barcode35.sorted.bam -b primerpair_pol.bed -p test -q 1 -s 4 > log.file

grep "Found" -A 10000 log.file | grep "primers in BED file" -A 100000 > ivar.log
