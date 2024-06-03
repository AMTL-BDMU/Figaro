process fastqc {
    container 'ufuomababatunde/nanoq:v0.10.0'

    tag "Check quality of ${sample}"

    publishDir path: { type == 'raw' ? "${params.outputDir}/01_fastqcRaw" : "${params.outputDir}/02_fastqcTrimmed" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fastq)
    val(type)

    output:
    tuple val(sample), path("*fastqc*"), emit: fastqc

    script:
    """
    fastqc --outdir . ${fastq}
    """
}
