// Creates a pdf report of the drug resistance mutation

process reportDrugResistance{
	container 'ufuomababatunde/rmarkdown:1.0.0'

	tag "Doing magic on $sample"

	
	publishDir (
	path: "${params.out_dir}/04_reports/",
	mode: 'copy',
	overwrite: 'true'
	)


	input:
	tuple val(sample), path(json)
	path(reportPDF)


	output:
	tuple val(sample), path('*.pdf'), optional: true, emit: report

	//4 spaces should be used instead of tab 
	script:
	"""
        Rscript -e 'rmarkdown::render("${reportPDF}", 
            params=list(
                mutation_comments="${params.sierraMutationDBComments}",
                dr_report_hivdb="${json}",
                output_file="hivdr_${sample}.pdf", output_dir = getwd())'
		
	"""
}
