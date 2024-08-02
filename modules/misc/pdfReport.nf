process pdfReport{
        container 'ufuomababatunde/rmarkdown:1.0.0'

        tag "Doing magic on $sample"

        
        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(json)
        path(markdownFile)


        output:
        tuple val(sample), path('*.pdf'), emit: report

        script:
        """
        Rscript -e 'rmarkdown::render("${markdownFile}", 
            params=list(
                logo="${params.logo}",
                dr_report_hivdb="${json}"),
                output_file="hivdr_${sample}.pdf", output_dir = getwd())'
        """
}


process pdfReport2{
        container 'ufuomababatunde/rmarkdown:1.0.0'

        tag "Doing magic on $sample"

        
        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(json)
        path(markdownFile)


        output:
        tuple val(sample), path('*.pdf'), emit: report

        script:
        """
        Rscript -e 'rmarkdown::render("${markdownFile}", 
            params=list(
                header="${params.reportHeader}",
                dr_report_hivdb="${json}"),
                output_file="hivdr_${sample}.pdf", output_dir = getwd())'
        """
}