process mosdepth {
	container 'nanozoo/mosdepth:0.3.2--892ca95'

	tag "checking coverage of $sample"

	publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}",
	mode: "copy",
	overwrite: "true"
	)

	
	input:
	tuple val(sample), path(bam), path(bai)

	output:
	path "*.dist.txt", emit: distribution
	tuple val(sample), path("*.txt"), emit: summary
    tuple val(sample), path("*.regions.bed"), emit: perBaseDepth

	script:
	"""
	mosdepth -b 1 -x ${sample} ${bam}

    gunzip ${sample}.regions.bed.gz
	"""
}