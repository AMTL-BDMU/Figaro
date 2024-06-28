process minimap2 {
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
            -O $params.minimapOpenPenalty \\
            -E $params.minimapExtenPenalty \\
            -a \\
            -t $params.thread \\
            $reference \\
            $fastq \\
            > ${sample}.sam
        """
}


process minimap2Iterate {
        container 'staphb/minimap2:2.28'

        tag "Assembling ${sample}"


        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq), path(reference)

        output:
        tuple val(sample), path("*sam"), emit: sam


        script:
        """
        minimap2 \\
            -O $params.minimapOpenPenalty \\
            -E $params.minimapExtenPenalty \\
            -a \\
            -t $params.thread \\
            $reference \\
            $fastq \\
            > ${sample}.sam
        """
}