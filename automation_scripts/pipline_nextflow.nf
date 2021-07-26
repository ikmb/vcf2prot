#!/usr/bin/env nextflow
//---------------------------------
// Author: Hesham ElAbd
// Copyright: Institute of Clinical Molecular Biology, 2021, Kiel, Germany. 
// Brief: A NextFlow script for parallizing the execution across multiple compute node, where different files will be processed on parallel
//----------------------------------
// introducing the script parameters
//----------------------------------
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Brief: 
    A wapper for the python automation script, can be used to parallize the execution of multiple file using a cluster of computional nodes
    Usage:
    The typical command for running the pipeline is as follows:
    Arguments:
      --input_dir [dir]                         The path to a directory contain a collection of bcf or bcf.gz files for execution 
      --reference [file]                        The path to a reference fasta file used for execution, i.e. generating personalized proteomes
      --output_dir [file]                       The directroy where results will be written, a directory for each input_file to hold the bcf or bcf.gz results there 
      --temp_dir [file]                         The path to create a temp directory where temp results will be created
      --batch_size [int]                        Number of probands to be processed in one batch of data
      --workers_per_node [int]                  Number of processes used by each node to compute results on parallel 
      -profile [str]                            Configuration profile to use. Can use multiple (comma separated)
                                                Available: docker, singularity, test, awsbatch and more
    """.stripIndent()
}
params.input_dir
params.reference
params.output_dir
params.temp_dir
params.batch_size
params.workers_per_node
// check that input is valid 
//--------------------------
params.input_dir = params.input_dir ?: {log.error "The path to the input directory has not been provided"; exit 1}()
params.reference = params.reference ?: {log error "The path to the reference has not been provided"; exit 1}()
params.output_dir = params.output_dir ?: {log.warn "The path to the write results has not been provided, results will be written to ./results"; return "./results"}()
params.temp_dir = params.temp_dir ?: {log.warn "The path to temp directory has not been provided, temp files will be created at ./temp"; return "./temp"}()
params.batch_size = params.batch_size ?: {log.warn "Batch size has not been provided, a default value of 1024 will be used"; return 1024}()
params.workers_per_node = params.workers_per_node? : {log.warn "Number of workers per node has not been provided, a default value of 16 will be used"; return 16}()
// Define the input channels 
//--------------------------
files_channel=Channel.fromPath(params.input_dir+'/*bcf*')
//--------------------------
// Define the process of the pipelines 
//------------------------------------
proccess RunPPGG{
    publishDir params.output_dir
    input:
    file bcf from files_channel 
    val ref from params.reference
    val outdir from params.output_dir
    val temp_dir from params.temp_dir
    val batch_size from params.batch_size
    val worker_per_node from params.workers_per_node
    """
    ./parallization_python -f ${bcf} -r ${ref} -o ${outdir} -t ${temp_dir} -b ${batch_size} -w ${worker_per_node} -k
    """
}
