process fastqc {
        container 'ufuomababatunde/nanoq:v0.10.0'

        tag "Check quality of ${sample}"


        publishDir {
            if (type == 'raw') {
                path: "${params.outputDir}/01_fastqcRaw"
            } else if (type == 'trimmed') {
                path: "${params.outputDir}/02_fastqcTrimmed"
            }
            mode: 'copy'
            overwrite: 'true'
        }


        input:
        tuple val(sample), path(fastq)
        value(type)

        output:
        tuple val(sample), path("*fastqc*"), emit: fastqc


        script:
        """
        fastqc \
        --outdir .
        $fastq
        """
}
