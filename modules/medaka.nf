process medakaPrelim {
        container 'ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/06_medakaPrelim",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(bam), path(bai)

        output:
        tuple val(sample), path("*.consensus.fasta"), emit:consensus

        script:
        """
        medaka consensus \\
            ${bam} \\
            ${sample}.hdf \\
            --model $params.medakaModel


        medaka stitch \\
            ${sample}.hdf \\
            $params.reference \\
            ${sample}.consensus.fasta

        sed -i "/^>/ s/.*/>${sample}/" ${sample}.consensus.fasta
        """
}



process medakaFinal {
        container 'ufuomababatunde/medaka:v1.11.4'

        tag "Creating consensus: ${sample}"


        publishDir (
        path: "${params.outDir}/08_medaka",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fasta), path(fastq)

        output:
        tuple val(sample), path("*.consensus.fasta"), emit: consensus


        script:
        """
        medaka_consensus \
            -t $params.thread \
            -m $params.medakaModel \
            -i $fastq \
            -d $fasta \
            -o medaka_dir

        mv medaka_dir/consensus.fasta ${sample}.consensus.fasta
        """
}
