process concatenate {
        tag "concatenating fastq.gz files in  ${samplePath}"


        publishDir (
        path: "${params.outputDir}/01_concatenated",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(samplePath)

        output:
        tuple val(sample), path ("*.fastq.gz"), emit: cat_fastq


        script:
        """
        cat ${samplePath}/*.fastq.gz > ${sample}.fastq.gz
        """
}