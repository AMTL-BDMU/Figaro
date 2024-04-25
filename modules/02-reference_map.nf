process minimap2 {
        cpus 1
        container 'staphb/minimap2'

        tag "mapping ${sample} to HBX2 POL reference"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(filtered_fastq), path (barcodePath)
        path(reference)

        output:
        tuple val(sample), path ("*.aligned.sam"), emit: ref_aligned_sam

        script:
        """
        minimap2 -t 8 -x map-ont -O 4 -E 2 -a ${reference} ${filtered_fastq} > ${sample}.aligned.sam
        """
}


process samtools {
        cpus 1
        container 'sanogenetics/htslib-bcftools-samtools'

        tag "samtools on ${sample}"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(sam)

        output:
        tuple val(sample), path ("*.sorted.bam"), path ("*.sorted.bam.bai"), emit: ref_aligned_bam

        script:
        """
        samtools view -bS ${sam} > ${sample}.aligned.bam
        samtools sort -o ${sample}.sorted.bam ${sample}.aligned.bam 
        samtools index -b ${sample}.sorted.bam
        """
}


process nanopolish {
        cpus 2
        container 'nanozoo/nanopolish:0.13.2--9069b5c'

        tag "making consensus sequence of ${sample}"


        publishDir (
        path: "${params.out_dir}/02_consensus",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(filtered_fastq), path (barcodePath), path(bam), path(bam_bai)
        path(reference)

        output:
        tuple val(sample), path ("*.fastq.index"), path ("*.fastq.index.fai"), path ("*.fastq.index.gzi"), path ("*.fastq.index.readdb"), path('*consensus.fasta'), emit: first_consensus

        script:
        """
        nanopolish index -d ${barcodePath}/../../fast5_pass/${sample} ${filtered_fastq}

        nanopolish variants --max-haplotypes 3000 --faster --fix-homopolymers --snps -t 32 -p 1 --reads ${filtered_fastq} --bam ${bam} --genome ${reference} > ${sample}.vcf

        nanopolish vcf2fasta -g ${reference} ${sample}.vcf > ${sample}.first_consensus.fasta
        """
}
