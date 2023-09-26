#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = "/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_convert/mzML"
params.input_datasets_folder = "/home/user/research/MSGFPLUS_search_workflow/data/search_td.fasta"


TOOL_FOLDER = "$baseDir/bin"


process MSGFPLUS_search {

    //container 'quay.io/biocontainers/thermorawfileparser:1.2.0--0'
    maxForks 10

    conda "$TOOL_FOLDER/conda_env.yml"

    memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 4
    publishDir "./nf_output", mode: 'copy'

    input:
    file spectra_file
    path dataset_file

    output:
    file '*.mzid'

    script:
    """
    java -Xmx3500M -jar $TOOL_FOLDER/MSGFPlus/MSGFPlus.jar -s ${spectra_file} -d ${dataset_file} -decoy XXX -o ${spectra_file.baseName}.mzid -t 30ppm
    """
}




workflow {
    spectra_files_ch = Channel.fromPath(params.input_spectra_folder+"/*.mzML")
    dataset_ch = params.input_datasets_folder
    idMZ_ch = MSGFPLUS_search(spectra_files_ch,dataset_ch)
}