process minimap2 {
        container 'staphb/minimap2:2.28'

        tag "Assembling ${sample}"


        publishDir (
        path: "${params.outDir}/03_minimap2",
        mode: 'copy',
        overwrite: 'true'
        )

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