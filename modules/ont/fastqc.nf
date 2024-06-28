process fastqcRaw {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of ${sample}"

        publishDir (
        path: "${params.outDir}/${workflow.name}/${name}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualRaw


        script:
        """
        fastqc --outDir . $fastq
        """
}


process fastqcTrimmed {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of ${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualTrimmed


        script:
        """
        fastqc --outDir . $fastq
        """
}