process consens_bcf {
	container 'biocontainers/python3-fast5:v0.5.8-1b2-deb_cv1'

	tag "Creating consensus of ${sample}"

	publishDir (
	path: "${params.out_dir}/03_consensus_bcf/",
	mode: 'copy',
	overwrite: 'true'
	)


	input:
	tuple val(sample), path(consensus), path(vcf), path(bgzipvcf), path(vcftbi)

	output:
	tuple val(sample), path("*.bcf.mixedSites.fasta"), emit: consens_bcf

	script:
	"""
	parse_vcf_bcf.py \
	--vcf ${bgzipvcf} \
	--out ${sample}.variants.tsv


	replaceNuc.py \
	--tsv ${sample}.variants.tsv \
	--infasta ${consensus} \
	--outfasta ${sample}.bcf.mixedSites.fasta
	"""
}



process consens_quasi {
	container 'codeyourinfra/python3:bionic'

	tag "Creating consensus of ${sample}"

	publishDir (
	path: "${params.out_dir}/03_consensus_quasi/",
	mode: 'copy',
	overwrite: 'true'
	)


	input:
	tuple val(sample), path(consensus), path(vcf), path(bgzipvcf), path(vcftbi)


	output:
	tuple val(sample), path("*.quasi.mixedSites.fasta"), emit: consens_quasi

	script:
	"""
	parse_vcf_quasi.py \
	--vcf ${bgzipvcf} \
	--out ${sample}.variants.tsv


	replaceNuc.py \
	--tsv ${sample}.variants.tsv \
	--infasta ${consensus} \
	--outfasta ${sample}.quasi.mixedSites.fasta
	"""
}



process consens_varscan {
	container 'biocontainers/bcftools:v1.9-1-deb_cv1'

	tag "Creating consensus of ${sample}"

	publishDir (
	path: "${params.out_dir}/03_consensus_varscan/",
	mode: 'copy',
	overwrite: 'true'
	)


	input:
	tuple val(sample), path(consensus), path(vcf), path(bgzipvcf), path(vcftbi)

	output:
	tuple val(sample), path("*.varscan.mixedSites.fasta"), emit: consens_varscan

	script:
	"""
	cat ${consensus} | \
	bcftools consensus \
	--iupac-codes \
	${bgzipvcf} > ${sample}.varscan.mixedSites.fasta
	"""
}


process consens_lofreq {
	container 'biocontainers/bcftools:v1.5_cv1'

	tag "Creating consensus of ${sample}"

	publishDir (
	path: "${params.out_dir}/03_consensus_lofreq/",
	mode: 'copy',
	overwrite: 'true'
	)


	input:
	tuple val(sample), path(consensus), path(vcf), path(bgzipvcf), path(vcftbi)

	output:
	tuple val(sample), path("*.lofreq.mixedSites.fasta"), emit: consens_lofreq

	script:
	"""
	cat ${consensus} | \
	bcftools consensus \
	--iupac-codes \
	${bgzipvcf} > ${sample}.lofreq.mixedSites.fasta
	"""
}
