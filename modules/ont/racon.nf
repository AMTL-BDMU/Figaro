process racon {
        container 'staphb/racon:1.5.0'

        tag "Polishing ${sample}"


        publishDir (
        path: "${params.outdir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fasta), path(fastq), path(sam)

        output:
        tuple val(sample), path("*.racon.fasta"), emit: raconFasta


        script:
        """
        racon \\
            --match $params.raconMatchScore \\
            --mismatch $params.raconMismatchScore \\
            --gap $params.raconGapPenalty \\
            --window-length $params.raconWindowLen \\
            --threads $params.thread \\
            $fastq \\
            $sam \\
            $fasta \\
            > ${sample}.racon.fasta
        """
}
