// enable dsl2
nextflow.enable.dsl=2


// Define a help message
def helpMessage = """
Usage:
    nextflow run Figaro {--ontAmplicon, --illuminaShotgun} --inDir reads/ --outDir results

Options:
    --help                  Show this help message and exit
    --inDir                 Directory containing the raw reads
    --outDir                Directory of the results

Choose only one workflow.
    --ontAmplicon           Run using the ONT amplicon workflow
    --illuminaShotgun       Run using the Illumina shotgun workflow
"""


// Check if help parameter is invoked and display help message
if (params.help) {
    println helpMessage
    exit 0
}


// import subworkflows
include {ontAmplicon} from './workflows/ont-wf.nf'
include {illuminaShotgun} from './workflows/illumina-wf.nf'

illuminaShotgun

workflow {
    main:
        if (params.ontAmplicon) {
            ontAmplicon()
        }

        else if (params.illuminaShotgun) {
            illuminaShotgun()
        }

        else {
            println helpMessage
        }
}