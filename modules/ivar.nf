process ivarPrelim {
        container 'staphb/ivar:1.4.2'

        tag "Trimming primers out of ${sample}"

        publishDir (
        path: "${params.outputDir}/04_ivarTrimPrelim",
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

process ivarFinal {
        container 'staphb/ivar:1.4.2'

        tag "Trimming primers out of ${sample}"

        publishDir (
        path: "${params.outputDir}/08_ivarTrimFinal",
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