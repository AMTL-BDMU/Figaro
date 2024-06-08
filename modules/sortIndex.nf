process sortIndex {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Sorting and indexing ${sample}"

        publishDir (
        path: "${params.outputDir}/05_sortIndex",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(bam)

        output:
        tuple val(sample), path("*.primerTrimmed.sorted.bam*"), path("*.primerTrimmed.sorted.bam.bai*"), emit: bamBai

        script:
        """
        samtools sort \
        ${bam} \
        -o ${sample}.primerTrimmed.sorted.bam

        samtools index ${sample}.primerTrimmed.sorted.bam
        """
}