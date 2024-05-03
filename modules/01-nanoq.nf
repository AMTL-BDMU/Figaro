#!/usr/bin/env nextflow

process cat {
        tag "concatenating fastq.gz files in  ${barcodePath.name}"


        publishDir (
        path: "${params.out_dir}/01_preprocess/",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path barcodePath

        output:
        tuple val(barcodePath.name), path ("*.fastq.gz"), path (barcodePath), emit: cat_fastq


        script:
        """
        cat ${barcodePath}/*.fastq.gz > ${barcodePath.name}.fastq.gz
        """
}



process nanoq {
        cpus 1
        container 'smobed/nanoq:latest'
        tag "trimming and filtering quality ${sample}"


        publishDir (
        path: "${params.out_dir}/01_preprocess",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq), path (barcodePath)

        output:
        tuple val(sample), path("*.nanoq.fastq.gz"), path (barcodePath), emit: nanoq_fastq
        path "${sample}_nanoq_report.txt"

        script:
        """
        nanoq \
        --input ${fastq} \
        --min-len 200 \
        --min-qual 12 \
        --output ${sample}.nanoq.fastq.gz \
        --report ${sample}_nanoq_report.txt
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
        fastqc --noextract --nogroup -o 
        """
}