#!/usr/bin/env nextflow

process cat {
        tag "concatenating fastq.gz files in  ${samplePath}"


        publishDir (
        path: "${params.out_dir}/01_preprocess/01_catFastq",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path samplePath

        output:
        path ("*.fastq.gz"), emit: cat_fastq


        script:
        """
        cat ${samplePath}/*.fastq.gz > ${samplePath}.fastq.gz
        """
}



process nanoq {
        cpus 1
        container 'smobed/nanoq:latest'
        tag "trimming and filtering quality ${sample}"


        publishDir (
        path: "${params.out_dir}/01_preprocess/01_nanoqFastq",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (cat_fastq)

        output:
        tuple val(sample), path("*.nanoq.fastq.gz"), emit: nanoq_fastq
        path "${sample}_nanoq_report.txt"

        script:
        """
        nanoq \
        --input ${cat_fastq} \
        --min-len 20 \
        --min-qual 12 \
        --output ${cat_fastq.baseName}.nanoq.fastq.gz \
        --report ${cat_fastq.baseName}_nanoq_report.txt
        """
}

process fastqc {
        cpus 1
        container 'staphb/fastqc:latest'
        tag "Quality report for ${sample}"


        publishDir (
        path: "${params.out_dir}/01_preprocess",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq), path (barcodePath)

        output:


        script:
        """
        
        """
}

process multiqc {
        cpus 1
        container 'staphb/fastqc:latest'
        tag "Quality report for ${sample}"


        publishDir (
        path: "${params.out_dir}/01_preprocess",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq), path (barcodePath)

        output:


        script:
        """
        
        """
}