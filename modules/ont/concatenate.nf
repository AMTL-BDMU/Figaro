process concatenate {
        tag "concatenating fastq.gz files in ${sample}"


        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(samplePath)

        output:
        tuple val(sample), path("*.fastq.gz"), emit: concatFastq


        script:
        """
        cat ${samplePath}/*.fastq.gz > ${sample}.fastq.gz
        """
}