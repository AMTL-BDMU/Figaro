process trimPrimer {
        container 'staphb/ivar:1.4.2'

        tag "Trimming primers out of ${sample}"

        publishDir (
        path: "${params.outDir}/05_ivarTrim",
        pattern: "*.primerTrimmed.bam",
        mode: 'copy',
        overwrite: 'true'
        )

        publishDir (
        path: "${params.outDir}/05_ivarLog",
        pattern: "*ivar.log",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(bam)

        output:
        tuple val(sample), path("*.primerTrimmed.bam"), emit:trimmedBam
        tuple val(sample), path("*ivar.log"), emit:ivarLog

        script:
        """
        ivar trim \
            -b $params.primer \
            -p ${sample}.primerTrimmed \
            -i ${bam} \
            -q 1 \
            -s 4 1> trim.log

        grep "Found" -A 10000 trim.log | grep "primers in BED file" -A 100000 > ${sample}.ivar.log
        """
}
