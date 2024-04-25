process minimap2_3 {
        cpus 1
        container 'staphb/minimap2'

        tag "mapping ${sample} to the second consensus fasta"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(second_consensus_fasta), path (filtered_fastq), path (barcodePath)

        output:
        tuple val(sample), path ("*.third.aligned.sam"), emit: con2_aligned_sam

        script:
        """
        minimap2 -t 8 -x map-ont -O 4 -E 2 -a ${second_consensus_fasta} ${filtered_fastq} > ${sample}.third.aligned.sam
        """
}


process samtools_3 {
        cpus 1
        container 'sanogenetics/htslib-bcftools-samtools'

        tag "samtools sort on ${sample}"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(sam)

        output:
        tuple val(sample), path ("*.sorted.bam"), path ("*.sorted.bam.bai"), emit: con2_aligned_bam


        script:
        """
        samtools view -bS ${sam} > ${sample}.third.aligned.bam
        samtools sort -o ${sample}.third.sorted.bam ${sample}.third.aligned.bam 
        samtools index -b ${sample}.third.sorted.bam  
        """
}

process nanopolish_3 {
        cpus 2
        container 'nanozoo/nanopolish:0.13.2--9069b5c'

        tag "making final consensus sequence of ${sample}"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(filtered_fastq), path (barcodePath), path (fastq_index), path (fastq_index_fai), path (fastq_index_gzi), path (fastq_index_readdb), path(consensus_fasta), path(bam), path(bam_bai), path (second_consensus_fasta)
        path(reference)

        output:
        tuple val(sample), path('*consensus.fasta'), path ('*.vcf'), emit: final_consensus

        script:
        """
        nanopolish variants --max-haplotypes 3000 --faster --fix-homopolymers --snps -t 32 -p 1 --reads ${filtered_fastq} --bam ${bam} --genome ${reference} > ${sample}.vcf
        nanopolish vcf2fasta -g ${reference} ${sample}.vcf > ${sample}.final_consensus.fasta
        """
}

process bgzip {
        cpus 1
        container 'sanogenetics/htslib-bcftools-samtools'

        tag 'bgzip ${sample}'

        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(final_consensus), path(vcf)

        output:
        tuple val(sample), path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: bgzip_vcf

        script:
        """
        bgzip ${vcf}
        tabix -p vcf ${sample}.vcf.gz
        bcftools index ${sample}.vcf.gz
        """
}
