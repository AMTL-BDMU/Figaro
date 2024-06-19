process sierra {
        container 'ufuomababatunde/sierralocal:1.4.0'

        //errorStrategy 'ignore'
        
        tag "Creating JSON file for $sample"

        publishDir (
        path: "${params.outDir}/11_sierra",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fasta)

        output:
        tuple val(sample), path('consensus*.json'), optional: true, emit: json
        path('*.drugResistanceScore.csv'), optional: true, emit: csv

        script:
        """
        if [ -s $fasta ]; then
            # The file is not-empty.
            sierralocal $fasta \
                -o tmp.json \
                -xml ${params.sierraXML}
            checkJSON.py --json tmp.json --sample ${sample}
        else
            echo "Skipping since there is no consensus sequence"
        fi


        if [ -s consensus_${sample}.json ]; then
        # If there is alignment.
            extractResistanceScore.py \
                --json consensus_${sample}.json \
                --csv ${sample}.drugResistanceScore.csv \
                --sample ${sample}
        else
            echo "Skipping since there is no enough alignment to any gene"
        fi
        """
}
