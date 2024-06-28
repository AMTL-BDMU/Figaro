process nanoq {
        container 'ufuomababatunde/nanoq:v0.10.0'

        tag "Filtering ${sample}"


        publishDir (
        path: "${params.outdir}/${task.process.replaceAll(":","_")}",
        pattern: "*.trimmed.fastq.gz",
        mode: 'copy',
        overwrite: 'true'
        )

        publishDir (
        path: "${params.outdir}/${task.process.replaceAll(":","_")}Report",
        pattern: "*.report.json",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*.trimmed.fastq.gz"), emit: trimmedFastq
        tuple val(sample), path("*.report.json"), emit: nanoqJSON


        script:
        """
        nanoq \\
            --input $fastq \\
            --json --report ${sample}.report.json \\
            --output ${sample}.trimmed.fastq.gz \\
            --min-len $params.nanoqMinReadLength \\
            --max-len $params.nanoqMaxReadLength \\
            --min-qual $params.nanoqMinReadQuality \\
            --max-qual $params.nanoqMaxReadQuality \\
            --trim-start $params.nanoqTrimStart \\
            --trim-end $params.nanoqTrimEnd
        """
}
