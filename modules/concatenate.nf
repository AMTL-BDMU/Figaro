process concatenate {
        tag "concatenating fastq.gz files in ${sample}"


        publishDir (
        path: "${params.outputDir}/01_concatenated",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(samplePath)

        output:
        tuple val(sample), path("*.fastq.gz"), emit: ch_concatFastq


        script:
        """
        cat ${samplePath}/*.fastq.gz > ${sample}.fastq.gz
        """
}