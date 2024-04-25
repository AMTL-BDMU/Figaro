#!/usr/bin/env nextflow

process cat {
        tag "concatenating fastq.gz files in  ${barcodePath.name}"


        publishDir (
        path: "${params.out_dir}/01_nanoq/",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path barcodePath

        output:
        path "${barcodePath.name}.fastq.gz", emit: cat_fastq


        script:
        """
        cat ${barcodePath}/*.fastq.gz > ${barcodePath.name}.fastq.gz
        """
}



process nanoq {
        cpus 1
//        container 'lanadelrea/nanoq:0.10.0'
//        container 'jimmyliu1326/nanoq'
        tag "trimming and filtering quality ${cat_fastq}"


        publishDir (
        path: "${params.out_dir}/01_nanoq",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(cat_fastq)

        output:
        path "${cat_fastq.SimpleName}.nanoq.fastq", emit: nanoq_fastq


        script:
        """
        nanoq --version
        nanoq \
        --input ${cat_fastq} \
        --min-len 20 \
        --min-qual 10 \
        --max-qual 20 \
        --trim-start 25 \
        --trim-end 25 \
        --output ${cat_fastq.SimpleName}.nanoq.fastq
        """
}
