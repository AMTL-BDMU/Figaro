process concatenate {
        tag "concatenating fastq.gz files in  ${samplePath}"


        publishDir (
        path: "${params.out_dir}/01_concatenated",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path(samplePath)

        output:
        path ("*.fastq.gz"), emit: cat_fastq


        script:
        """
        cat ${samplePath}/*.fastq.gz > ${samplePath}.fastq.gz
        """
}