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

        # Set overallResult default as "Failed"
        # This will be updated later on the pipeline
        update_json.py \\
            --json starting_blank.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature overallResult \\
            --value Failed

        update_json.py \\
            --json ${sample}.numlen.updated.json \\
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
            --feature meanReadLength_initial \\
            --value \${readMeanLength_raw}

        update_json.py \\
            --json ${sample}.numlen.updated.json \\
            --out ${sample}.numlen.updated.json \\
            --sample ${sample} \\
            --feature meanReadLength_final \\
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


process getMappedPercent {
        container 'ufuomababatunde/samtools:v1.17-pymodulesv4-seqtkv1.4-r130-dirty'

        tag "${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(sam), path(json)

        output:
        tuple val(sample), path("*.mapped.updated.json"), emit: mappedJSON

        script:
        """      
        mappedPercent=\$(samtools flagstat ${sam} -O tsv | grep "mapped %" | head -n1 | awk -F "\\t" '{print \$1}' | tr -d %)

        update_json.py \\
            --json ${json} \\
            --out ${sample}.mapped.updated.json \\
            --sample ${sample} \\
            --feature mappedReads \\
            --value \${mappedPercent}


        """
}





process getBAMinfo {
//        container 'ufuomababatunde/samtools:v1.17-pymodulesv4-seqtkv1.4-r130-dirty'
        conda '/apps/miniconda3/envs/figaro'
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

        meanCoverage=\$(samtools coverage ${bam} --no-header | awk -F "\\t" '{print \$6}')
        meanDepth=\$(samtools coverage ${bam} --no-header | awk -F "\\t" '{print \$7}')


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

        update_json.py \\
            --json ${sample}.depth.updated.json \\
            --out ${sample}.depth.updated.json \\
            --sample ${sample} \\
            --feature genomeCoverage \\
            --value \${meanCoverage}

        update_json.py \\
            --json ${sample}.depth.updated.json \\
            --out ${sample}.depth.updated.json \\
            --sample ${sample} \\
            --feature meanDepth \\
            --value \${meanDepth}
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
        val(sampleList)

        output:
        path("filledIn_combined.json"), emit: combinedJSON
        path("report.html"), emit: htmlReport

        script:
        """
        combine_json.py \\
            --input ${json} \\
            --output combined.json
        
        check_jsonFeatures.py \\
            --inJSON combined.json \\
            --outTXT missingFeatures.txt \\
            --outJSON filledIn_combined.json

        update_overallResult.py \\
            --inJSON filledIn_combined.json \\
            --outJSON updated_overAllResult.json \\
            --result Passed \\
            --sample ${sampleList.join(' ')}


        cp ${params.templateHtmlRunReport} .

        generate_report.py \\
            --json updated_overAllResult.json \\
            --template report_template.html \\
            --report report.html
        """
}