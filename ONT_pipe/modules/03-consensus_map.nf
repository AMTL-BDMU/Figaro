process minimap2_2 {
        cpus 1
        container 'staphb/minimap2'

        tag "mapping ${sample} to the first consensus fasta"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(first_consensus_fasta), path(filtered_fastq), path (barcodePath)

        output:
        tuple val(sample), path ("*.second.aligned.sam"), emit: con1_aligned_sam

        script:
        """
        minimap2 -t 8 -x map-ont -O 2 -E 1 -a ${first_consensus_fasta} ${filtered_fastq} > ${sample}.second.aligned.sam
        """
}


process samtools_2 {
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
        tuple val(sample), path ("*.sorted.bam"), path ("*.sorted.bam.bai"), emit: con1_aligned_bam


        script:
        """
        samtools view -bS ${sam} > ${sample}.second.aligned.bam
        samtools sort -o ${sample}.second.sorted.bam ${sample}.second.aligned.bam 
        samtools index -b ${sample}.second.sorted.bam  
        """
}

process nanopolish_2 {
        cpus 2
        container 'nanozoo/nanopolish:0.13.2--9069b5c'

        tag "making consensus sequence of ${sample}"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path (filtered_fastq), path (barcodePath), path (fastq_index), path (fastq_index_fai), path (fastq_index_gzi), path (fastq_index_readdb), path (consensus_fasta), path(bam), path(bam_bai)
        path(reference)

        output:
        tuple val(sample), path('*consensus.fasta'), emit: second_consensus

        script:
        """
        nanopolish variants --max-haplotypes 3000 --faster --fix-homopolymers --snps -t 32 -p 1 --reads ${filtered_fastq} --bam ${bam} --genome ${reference} > ${sample}.vcf

        nanopolish vcf2fasta -g ${reference} ${sample}.vcf > ${sample}.second_consensus.fasta
        """
}

