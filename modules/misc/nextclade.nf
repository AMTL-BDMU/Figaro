process nextcladeDB {
    container 'ufuomababatunde/nextclade:3.8.2'

    tag "Downloading neherlab/hiv-1 dataset"

    output:
    path "hiv1_db", emit: database

    script:
    """
    nextclade dataset get \\
        --name "neherlab/hiv-1" \\
        --output-dir hiv1_db
    """
}


process nextcladeSubtype {
        container 'ufuomababatunde/nextclade:3.8.2'

        errorStrategy 'ignore'
        
        tag "Subtyping $sample"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )


        input:
        tuple val(sample), path(fasta)
        each file(hiv1_db)

        output:
        tuple val(sample), path('*tsv'), emit: tsv

        script:
        """
        nextclade run \\
            --input-dataset ${hiv1_db} \\
            --output-tsv ${sample}.tsv \\
            ${fasta}

        
        """
}
