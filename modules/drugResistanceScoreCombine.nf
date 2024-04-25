// Combine csv files containing drug resistance scores from hydra 



process drugResistanceScoreCombine {
	container 'ufuomababatunde/biopython:v1.2.0'

	tag "Collecting samples"

	
	publishDir (
	path: "${params.out_dir}/04_sierra/",
	mode: 'copy',
	overwrite: 'true'
	)

	//errorStrategy 'ignore'
	
	input:
	path(csv)

	output:
	path("combinedDRscores.csv")


	script:
	"""
	combineDRscores.py \
	combinedDRscores.csv \
	$csv

	"""
}
