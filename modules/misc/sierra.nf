process sierra {
        container 'ufuomababatunde/sierralocal:0.4'

        errorStrategy 'ignore'
        
        tag "Creating JSON file for $sample"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fasta), path(tsv)

        output:
        tuple val(sample), path('subtype*.json'), optional: true, emit: json
        path('*.drugResistanceScore.csv'), optional: true, emit: csv

        script:
        """
        if [ -s $fasta ]; then
            # The file is not-empty.
            sierralocal $fasta \\
                -o tmp.json \\
                -xml ${params.sierraXML}
                
            checkJSON.py --json tmp.json --sample ${sample}
        else
            echo "Skipping since there is no consensus sequence"
        fi


        if [ -s consensus_${sample}.json ]; then
        # If there is alignment.
            extractResistanceScore.py \\
                --json consensus_${sample}.json \\
                --csv ${sample}.drugResistanceScore.csv \\
                --sample ${sample}
        else
            echo "Skipping since there is no enough alignment to any gene"
        fi


        # Add subtype from nextclade
        subtype=\$(awk -F "\\t" '{print \$3}' ${tsv} | tail -n1)

        addSubtype.py \\
            --inJSON consensus_${sample}.json \\
            --subtype \${subtype}
            --outJSON subtype_${sample}.json
        """
}
