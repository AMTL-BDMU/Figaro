process downloadRef {
        tag "downloading HBX2 POL reference"


        publishDir (
        path: "${params.out_dir}/02_referencemap",
        mode: 'copy',
        overwrite: 'true'
        )

        output:
        path "hxb2_pol.fas", emit: reference

        script:
        """
        wget https://raw.githubusercontent.com/phac-nml/quasitools/master/quasitools/data/hxb2_pol.fas
        """
}


process minimap2 {
        cpus 1
        container 'staphb/minimap2'

        tag "mapping ${nanoq_fastq} to HBX2 POL reference"


        publishDir (
        path: "${params.out_dir}/02_referencemap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(nanoq_fastq)
        file(reference)

        output:
        path "${nanoq_fastq.SimpleName}.aligned.sam", emit: ref_aligned_sam

        script:
        """
        minimap2 -t 8 -ax map-ont ${reference} ${nanoq_fastq} > ${nanoq_fastq.SimpleName}.aligned.sam
        """
}


process samtools {
        cpus 1
        container 'ontresearch/medaka'

        tag "making consensus from ${ref_aligned_sam}"


        publishDir (
        path: "${params.out_dir}/02_referencemap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(ref_aligned_sam)
        file(reference)

        output:
        path "${ref_aligned_sam.SimpleName}.sorted.bam", emit: ref_sorted_bam
        path "${ref_aligned_sam.SimpleName}.sorted.bam.bai"
        path "${ref_aligned_sam.SimpleName}.first.fasta", emit: first_fasta

        script:
        """
        samtools view -bS ${ref_aligned_sam} > ${ref_aligned_sam.SimpleName}.aligned.bam
        samtools sort -o ${ref_aligned_sam.SimpleName}.sorted.bam ${ref_aligned_sam.SimpleName}.aligned.bam 
        samtools index -b ${ref_aligned_sam.SimpleName}.sorted.bam
        bcftools mpileup -f ${reference} ${ref_aligned_sam.SimpleName}.sorted.bam > ${ref_aligned_sam.SimpleName}.mpileup
        bcftools call -mv -O b -o ${ref_aligned_sam.SimpleName}.bcf ${ref_aligned_sam.SimpleName}.mpileup
        bcftools view -O v -o ${ref_aligned_sam.SimpleName}.vcf ${ref_aligned_sam.SimpleName}.bcf
        bgzip ${ref_aligned_sam.SimpleName}.vcf
        tabix -p vcf ${ref_aligned_sam.SimpleName}.vcf.gz
        bcftools consensus -f ${reference} ${ref_aligned_sam.SimpleName}.vcf.gz > ${ref_aligned_sam.SimpleName}.first.fasta 
        """
}


process medaka {
        cpus 2
        container 'ontresearch/medaka'

        tag "polishing consensus sequence of ${first_fasta}"


        publishDir (
        path: "${params.out_dir}/02_referencemap",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        file(first_fasta)
        file(reference)

        output:
        tuple val(first_fasta.SimpleName), path('*consensus.fasta'), emit: first_fasta_polished

        script:
        """
        medaka_consensus -i ${first_fasta} -d ${reference} -o ${first_fasta.SimpleName}_medaka -t 8 -m r941_min_hac_g507

        cp ${first_fasta.SimpleName}_medaka/consensus.fasta ${first_fasta.SimpleName}.consensus.fasta
        """
}
