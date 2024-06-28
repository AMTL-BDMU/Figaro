process medaka {
        container 'ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/06_medaka",
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
