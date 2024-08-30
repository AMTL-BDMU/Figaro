process hydra {
	container 'ufuomababatunde/quasitools:v0.7.0'

	tag "Doing magic on $sample"

	
	publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}",
    pattern: "*fasta",
	mode: 'copy',
	overwrite: 'true'
	)


    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}Alignment",
    pattern: "*.align.bam*",
    mode: 'copy',
    overwrite: 'true'
    )


	input:
	tuple val(sample), path(fastq_1), path(fastq_2)

	output:
	tuple val(sample), path("*.consensus.fasta"), emit: consensus
    tuple val(sample), path("*.align.bam"), path("*.align.bam.bai"), emit: bamBai


	script:
	"""
	quasitools hydra ${fastq_1} ${fastq_2} \\
        --generate_consensus \\
        --mask_reads \\
        --reporting_threshold ${params.hydraReportThreshold} \\
        --min_read_qual ${params.hydraMinReadQuality} \\
        --consensus_pct ${params.hydraConsensusPercent} \\
        --length_cutoff ${params.hydraMinReadLength} \\
        --score_cutoff ${params.hydraMinScoreCutoff} \\
        --min_variant_qual ${params.hydraMinVariantQuality} \\
        --min_dp ${params.hydraMinVariantDepth} \\
        --min_ac ${params.hydraMinAlleleCount} \\
        --min_freq ${params.hydraMinVariantFrequency} \\
        --id ${sample} \\
        --output_dir ${sample}

    cp ${sample}/consensus.fasta ${sample}.consensus.fasta

    cp ${sample}/align.bam ${sample}.align.bam
    cp ${sample}/align.bam.bai ${sample}.align.bam.bai
	"""

}