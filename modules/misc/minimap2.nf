process minimap2SE {
        container 'staphb/minimap2:2.28'

        tag "Assembling ${sample}"


        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)
        path(reference)

        output:
        tuple val(sample), path("*sam"), emit: sam


        script:
        """
        minimap2 \\
            -ax $params.minimapPreset \\
            -t $params.thread \\
            $reference \\
            $fastq \\
            > ${sample}.sam
        """
}


process minimap2PE {
        container 'staphb/minimap2:2.28'

        tag "Assembling ${sample}"


        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq_1), path(fastq_2)
        path(reference)

        output:
        tuple val(sample), path("*sam"), emit: sam


        script:
        """
        minimap2 \\
            -ax $params.minimapPreset \\
            -t $params.thread \\
            $reference \\
            $fastq_1  $fastq_2 \\
            > ${sample}.sam
        """
}