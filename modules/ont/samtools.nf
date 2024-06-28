process sam2bam {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Converting sam to bam: ${sample}"

        publishDir (
        path: "${params.outDir}/04_sam2bam",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(sam)

        output:
        tuple val(sample), path("*.bam"), emit: bam

        script:
        """
        samtools view \
            -bS ${sam} > ${sample}.bam
        """
}


process sortIndex {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Sorting and indexing ${sample}"

        publishDir (
        path: "${params.outDir}/04_sortIndex",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(bam)

        output:
        tuple val(sample), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bamBai

        script:
        """
        samtools sort \
            ${bam} \
            -o ${sample}.sorted.bam

        samtools index ${sample}.sorted.bam
        """
}
