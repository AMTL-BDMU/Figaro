process minimapPrelim {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Aligning ${sample}"

        publishDir (
        path: "${params.outputDir}/03_minimapPrelim",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)
        path(reference)

        output:
        tuple val(sample), path("*.sorted.bam"), emit: bam

        script:
        """
        minimap2 \
            -a ${reference} \
            ${fastq} | samtools sort - -o ${sample}.sorted.bam
        """
}

process minimapFinal {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Aligning ${sample}"

        publishDir (
        path: "${params.outputDir}/07_minimapFinal",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)
        path(reference)

        output:
        tuple val(sample), path("*.sorted.bam"), emit: bam

        script:
        """
        minimap2 \
            -a ${reference} \
            ${fastq} | samtools sort - -o ${sample}.sorted.bam
        """
}