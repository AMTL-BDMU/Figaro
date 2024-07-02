    process fastP {
    container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'

    tag "trimming $sample"


    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("*1.fastq"), path("*2.fastq"), emit: trimmed
    tuple val(sample), path("*.json"), emit: fastP_json
    tuple val(sample), path("*.html"), emit: fastP_html


    script:
    """
    fastp \\
        -q $params.fastPPhred \\
        -i $fastq_1 \\
        -I $fastq_2 \\
        -o ${sample}.trimmed_R1.fastq \\
        -O ${sample}.trimmed_R2.fastq \\
        -j ${sample}.fastp.json \\
        -h ${sample}.fastp.html
    """

    }