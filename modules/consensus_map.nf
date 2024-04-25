process minimap2_2 {
        cpus 1
        container 'staphb/minimap2'

        tag "mapping ${nanoq_fastq} to the first consensus fasta"


        publishDir (
        path: "${params.out_dir}/03_consensusmap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fasta)
        file(nanoq_fastq)

        output:
        path "${sample}.second.aligned.sam", emit: con1_aligned_sam

        script:
        """
        minimap2 -t 8 -ax map-ont ${fasta} ${nanoq_fastq} > ${sample}.second.aligned.sam
        """
}


process samtools_2 {
        cpus 1
        container 'ontresearch/medaka'

        tag "samtools sort on ${con1_aligned_sam}"


        publishDir (
        path: "${params.out_dir}/03_consensusmap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(con1_aligned_sam)
        tuple val(sample), path(fasta)

        output:
        path "${sample}.second.sorted.bam"
        path "${sample}.second.sorted.bam.bai"
        path "${sample}.second.fasta", emit: second_fasta

        script:
        """
        samtools view -bS ${con1_aligned_sam} > ${sample}.second.aligned.bam
        samtools sort -o ${sample}.second.sorted.bam ${sample}.second.aligned.bam 
        samtools index -b ${sample}.second.sorted.bam
        bcftools mpileup -f ${fasta} ${sample}.second.sorted.bam > ${sample}.second.mpileup
        bcftools call -mv -O b -o ${sample}.second.bcf ${sample}.second.mpileup
        bcftools view -O v -o ${sample}.second.vcf ${sample}.second.bcf
        bgzip ${sample}.second.vcf
        tabix -p vcf ${sample}.second.vcf.gz
        bcftools consensus -f ${fasta} ${sample}.second.vcf.gz > ${sample}.second.fasta 
        """
}


process medaka_2 {
        cpus 2
        container 'ontresearch/medaka'

        tag "polishing consensus sequence of ${second_fasta}"


        publishDir (
        path: "${params.out_dir}/03_consensusmap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(second_fasta)
        tuple val(sample), path(first_fasta)

        output:
        tuple val(second_fasta.SimpleName), path('*.second.consensus.fasta'), emit: second_fasta_polished

        script:
        """
        medaka_consensus -i ${second_fasta} -d ${first_fasta} -o ${second_fasta.SimpleName}_medaka_second -t 8 -m r941_min_hac_g507

        cp ${second_fasta.SimpleName}_medaka_second/consensus.fasta ${second_fasta.SimpleName}.second.consensus.fasta
        """
}
