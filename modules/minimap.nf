process minimapPrelim {
        container 'ufuomababatunde/minimap2:v2.26-samtoolsv1.18'

        tag "Aligning ${sample}"

        publishDir (
        path: "${params.outDir}/03_minimapPrelim",
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
            -O $params.minimapOpenPenalty \
            -E $params.minimapExtenPenalty \
            -a ${reference} \
            ${fastq} | samtools sort - -o ${sample}.sorted.bam
        """
}


process minimapFinal {
        container 'staphb/minimap2:2.28'

        tag "Assembling ${sample}"


        // publishDir (
        // path: "${params.outDir}/04_minimap2",
        // mode: 'copy',
        // overwrite: 'true'
        // )

        input:
        tuple val(sample), path(fasta), path(fastq)

        output:
        tuple val(sample), path("*sam"), emit: sam


        script:
        """
        minimap2 \
            -O $params.minimapOpenPenalty \
            -E $params.minimapExtenPenalty \
            -a \
            -t $params.thread \
            $fasta \
            $fastq \
            > ${sample}.sam
        """
}
