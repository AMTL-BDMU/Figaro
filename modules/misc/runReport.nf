process getReadNumberLength {
        container 'ufuomababatunde/seqkit-pymodule:v2.8.2'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq_raw), path(fastq_trimmed)

        output:
        tuple val(sample), path("*.numlen.updated.json"), emit: numlenJSON

        script:
        """
        readNumber_raw=\$(seqkit stats $fastq_raw | awk -F " " '{print \$4}' | tail -n1 | sed 's/,//g')
        readMeanLength_raw=\$(seqkit stats $fastq_raw | awk -F " " '{print \$7}' | tail -n1 | sed 's/,//g')

        readNumber_trimmed=\$(seqkit stats $fastq_trimmed | awk -F " " '{print \$4}' | tail -n1 | sed 's/,//g')
        readMeanLength_trimmed=\$(seqkit stats $fastq_trimmed | awk -F " " '{print \$7}' | tail -n1 | sed 's/,//g')

        passedReadsProp=\$(echo "\$readNumber_trimmed/\$readNumber_raw * 100" | bc -l)





        create_blankJSON.py \\
            --title ${params.outDir} \\
            --outJSON starting_blank.json


        update_json.py \\
            --json starting_blank.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature readsNumber_initial \\
            --value \${readNumber_raw}

        update_json.py \\
            --json ${sample}.numlen.updated.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature readsProportion_passed \\
            --value \${passedReadsProp}

        update_json.py \\
            --json ${sample}.numlen.updated.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature medianReadLength_initial \\
            --value \${readMeanLength_raw}

        update_json.py \\
            --json ${sample}.numlen.updated.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature medianReadLength_final \\
            --value \${readMeanLength_trimmed}
        """
}


process getPhredTrimmed {
        container 'ufuomababatunde/seqkit-pymodule:v2.8.2'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fastq), path(json)

        output:
        tuple val(sample), path("*.phred.updated.json"), emit: phredJSON

        script:
        """
        median_quality.py \\
            --inFastq ${fastq} \\
            --outSummary ${sample}.summaryQual.tsv
        
        readBP=\$(tail -n+2 ${sample}.summaryQual.tsv | awk -F "\\t" '{print \$1}' | paste -sd ',')

        readPhred=\$(tail -n+2 ${sample}.summaryQual.tsv | awk -F "\\t" '{print \$2}' | paste -sd ',')


        update_json.py \\
            --json ${json} \\
            --out ${sample}.phred.updated.json \\
            --sample ${sample} \\
            --feature qc_bp \\
            --value \${readBP}

        update_json.py \\
            --json ${sample}.phred.updated.json \\
            --out ${sample}.phred.updated.json \\
            --sample ${sample} \\
            --feature qc_phred \\
            --value \${readPhred}
        """
}


process getDepth {
        container 'ufuomababatunde/samtools:v1.17-pymodulesv4-seqtkv1.4-r130-dirty'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(bam), path(bai), path(json)

        output:
        path("*.depth.updated.json"), emit: depthJSON

        script:
        """      
        genomeBP=\$(samtools depth -a ${bam} | awk -F "\\t" '{print \$2}' | paste -sd ',')

        genomeDepth=\$(samtools depth -a ${bam} | awk -F "\\t" '{print \$3}' | paste -sd ',')


        update_json.py \\
            --json ${json} \\
            --out ${sample}.depth.updated.json \\
            --sample ${sample} \\
            --feature alignment_bp \\
            --value \${genomeBP}

        update_json.py \\
            --json ${sample}.depth.updated.json \\
            --out ${sample}.depth.updated.json \\
            --sample ${sample} \\
            --feature alignment_depth \\
            --value \${genomeDepth}
        """
}


process checkCombineJSON {
        container 'ufuomababatunde/seqkit-pymodule:v2.8.2'

        tag "Generating Run report"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path(json)

        output:
        path("filledIn_combined.json"), emit: combinedJSON

        script:
        """
        combine_json.py \\
            --input ${json} \\
            --output combined.json
        
        check_jsonFeatures.py \\
            --inJSON combined.json \\
            --outTXT missingFeatures.txt \\
            --outJSON filledIn_combined.json
        """
}



process htmlRunReport {
        container 'ufuomababatunde/seqkit-pymodule:v2.8.2'

        tag "Generating Run report"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        path(json)

        output:
        path("report.html"), emit: htmlReport

        script:
        """
        generate_report.py \\
            --json ${json} \\
            --template ${params.templateHtmlRunReport} \\
            --report report.html
        """
}