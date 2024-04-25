process codfreq {

        container 'hivdb/codfreq-runner:latest'

        tag "fastq to codfreq ${nanoq_fastq}"


        publishDir (
        path: "${params.out_dir}/01_codfreq",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(con1_aligned_sam)
        file(nanoq_fastq)

        output:
        path "", emit: 

        script:
        """
        fastq2codfreq         
        """
}
