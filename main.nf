// enable dsl2
nextflow.enable.dsl=2


// Define a help message
def helpMessage = """
Usage:
    nextflow run Figaro \
        --ont-amplicon \
        --inDir reads/ \
        --outDir results

Options:
    --help                  Show this help message and exit
    --ont-amplicon          Run using the ONT amplicon workflow
    --inDir                 Directory containing the raw reads
    --outDir                Directory of the results
"""


// Check if help parameter is invoked and display help message
if (params.help) {
    println helpMessage
    exit 0
}