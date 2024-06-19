process fastqcRaw {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of ${sample}"

        publishDir (
        path: "${params.outDir}/01_fastqcRaw",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualRaw


        script:
        """
        fastqc --outdir . $fastq
        """
}


process fastqcTrimmed {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of ${sample}"

        publishDir (
        path: "${params.outDir}/02_fastqcTrimmed",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualTrimmed


        script:
        """
        fastqc --outdir . $fastq
        """
}