#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: A Python script for extracting batches of patients from a bcf.gz files, write them to a VCF file, and then run them in parallel
using a pool of processes.
@copyright: Institute of clinical molecular biology (IKMB) 2021, Kiel, Germany. 
"""
## load the modules and install if not installed 
import argparse
import time
import os
import subprocess as sp 
from typing import List 
from concurrent import futures
import random 
try: 
    import pandas as pd 
except ModuleNotFoundError:
    os.system("pip install pandas")
import concurrent.futures
try: 
    from tqdm import tqdm 
except ModuleNotFoundError:
    os.system("pip install tqdm")
## record the start time
#-----------------------
start_time=time.ctime()
## define the user argument 
#--------------------------
parser=argparse.ArgumentParser(description="A wrapper for PPGG, can be used for extracting a batch of patient from a bcf.gz into a VCF file \
    and call PPGG on the generate batch to generated personalized proteome for each patient in the extracted batch. Different batches are executed concurrently\
        using a pool of processes")

parser.add_argument("-f,--file",help="The path to the input bcf.gz file", default='?')
parser.add_argument("-r,--ref", help="The path to the reference fasta file", default="?")
parser.add_argument("-o,--output_file", help="The path to write the generated fasta sequences, defaults to the current working directory",
                    defaults=os.getcwd())
parser.add_argument("-t,--temp_dir", help="The path to a temp directory to store VCF files, defaults to the current working directory",
                    defaults=os.getcwd())
parser.add_argument("-w,--num_workers", help="number of worker processes, defaults to number of available CPU cores", 
                            default=os.cpu_count())
parser.add_argument("-b, --batch_size",help="The number of patients in a batch, defaults to 100", default=100, type=int)
parser.add_argument("-e,--exe", help="The path to PPGG executables",default=None)
parser.add_argument('-k,--set_input_name_as_base',help="A boolean switch for controlling the writting output, if True, a directory with\
    same name as the input is created in the output directory and the temp directory, where files will be written.", default=False, action='store_true')
## Parse user provided arguments 
#-------------------------------
args=parser.parse_args()
# Check that the input file exist 
if args.file=="?":
    raise ValueError("Input file has not been provided")
elif not os.path.exists(args.file):
    raise ValueError(f"The provided path: {args.file} does not exist!!")
# Check that the reference file exist 
if args.ref=="?":
        raise ValueError("Reference Fasta file has not been provided")
elif not os.path.exists(args.file):
    raise ValueError(f"The provided path: {args.file} does not exist!!")
# the path of the executable to the path 
## define an analysis function 
def get_patient_name(path2file:str, path2temp:str)->None:
    """Write the name of patients in the provided bcf.gz file to a text file named 'name_probands.txt'
       in the temp directory

    Args:
        path2file (str): The path to the bcf.gz file containing all patients genomic alterations.
        path2temp (str): The path to the temp directory
    """
    try: 
        sp.run(f"bcftools query -l {path2file} > {os.path.join(path2temp,'name_probands.txt')}")
    except sp.SubprocessError as exp:
        raise ValueError(f"Getting the probands names failed with the following error code: {str(exp)}")
    return

def load_patient_name(path2file:str)->pd.DataFrame:
    """Load patient names into a pd.DataFrame

    Args:
        path2file (str): The path to a text file containing the name of each patient in the BCF compressed file

    Returns:
        pd.DataFrame: a table of the results
    """
    try: 
        names=pd.read_csv(path2file,header=None)
    except Exception as exp:
        raise IOError(f"Loading the patient names failed with the following error: {str(exp)}")
    return names

def extract_patient_from_bcf_compressed(path2file:str, patients: List[str], path2temp:str)->str:
    """ Extract a batch of samples from the compressed BCF file and write them on parallel to a BCF File 

    Args:
        path2file (str): the path the VCF file 
        patients (List[str]): a list containing the name of patient to extract their data from the file 
        path2temp (str): the path to the temp directory, where the vcf file will be written

    Returns:
        str: the absolute path to the generate VCF file 
    """
    # 1. write patient names to a text file in the temp directory 
    # 1.a generate a randome identifier 
    batch_identifer=str(random.randint(int(1e6),int(1e7)))
    # 1.b create a file name 
    file_path=os.path.join(path2temp,'patinet_batch_'+batch_identifer+'list_patient.txt')
    # 1.c write patient names on parallel
    try:
        with open(file_path,'w') as writer_file:
            for patient in patients:
                writer_file.write(patient+'\n')
    except Exception as exp: 
        raise IOError(f"Writing the patients file to the path: {file_path}, failed with the following error: {str(exp)}")
    # 2- extract a vcf file from a bcf file 
    result_file=os.path.join(path2temp,'patinet_batch_'+batch_identifer+'alterations.vcf')
    try:
        sp.run(f"bcftools view -c1 -Ov -S {os.path.abspath(file_path)} -o {os.path.abspath(result_file)} {os.path.abspath(path2file)}")
    except sp.SubprocessError as exp:
        raise RuntimeError(f"Extracting the batch of data, failed with the following error: {str(exp)}")
    return result_file

def worker_function(path2chunK:str, path2ref:str, result_path:str)->None:
    """Runs PPGG with the provided batch of probands on a separate process

    Args:
        path2chunK (str): The path to a vcf file containing the batch of samples used in the analysis 
        path2ref (str): The path to the reference FASTA file 
        result_path (str): The path to write the resulting FASTA file 

    Raises:
        RuntimeError: incase calling PPGG failed
    """
    try:
        if args.exe==None:
            sp.run(f"ppgg_rust -f {path2chunK} -r {path2ref}, -o {result_path} -g st",
                                stdout=sp.DEVNULL,check=True)
        else:
            sp.run(f".{args.exe}/ppgg_rust -f {path2chunK} -r {path2ref}, -o {result_path} -g st",
                                stdout=sp.DEVNULL,check=True)
    except sp.SubprocessError as exp:
        raise RuntimeError(f"Calling PPGG with the following chunk failed: {path2chunK}")

def personalization_task(path2file:str,patients:List[str], path2temp:str, path2ref:str, path2res:str)->None:
    """ Extract the genomic alteration of a batch of patients from a VCF file, write the results to a VCF file\
         and call PPGG to generate the personalized Fasta sequences. 

    Args:
        path2file (str): the path the VCF file. 
        patients (List[str]): a list containing the name of patient to extract their data from the file. 
        path2temp (str): The path to the temp directory, where the vcf file will be written.
        path2temp (str): The path to the reference FASTA file.
        path2res (str): The path to write the generated Fasta Files 
    """
    # 1. extract vcf files: 
    path2vcf=extract_patient_from_bcf_compressed(path2file,patients,path2temp)
    # 2. generating the results
    worker_function(path2vcf,path2ref,path2res)
    return
## Start the execution part of the script 
#----------------------------------------
# 1. Create a temp directory 
if os.path.exists(args.temp_dir) and not args.set_input_name_as_base:
    pass
else: 
    try:
        if args.set_input_name_as_base:
            base_dir=args.file.split('/')[-1].split('.')[0] 
            args.temp_dir=os.path.join(args.temp_dir,base_dir)
            os.mkdir(args.temp_dir)
        else:
            os.mkdir(args.temp_dir)
    except Exception as exp:
        raise RuntimeError(f"Creating a temp directory failed with the following error:{str(exp)}") 
# 2. Create an output directory 
if os.path.exists(args.output_dir) and not args.set_input_name_as_base:
    pass
else: 
    try:
        if args.set_input_name_as_base:
            base_dir=args.file.split('/')[-1].split('.')[0] 
            args.output_dir=os.path.join(args.output_dir,base_dir)
            os.mkdir(args.output_dir)
        else:
            os.mkdir(args.output_dir)
    except Exception as exp:
        raise RuntimeError(f"Creating a temp directory failed with the following error:{str(exp)}") 
# 3. Create a list of patient 
get_patient_name(args.input_file,args.temp_dir)
# 4. load the patient name as a list of names
patient_names=load_patient_name(os.path.join(args.temp_dir,'name_probands.txt')).iloc[:,0].to_list()
# 5. distribute the load among different processes 
# 5.1 create a list to hold future objects
jobs=[]
# 5.2 compute the load of each process
load=patient_names.shape[0]/args.batch_size
# 5.3 create a workers pool and submit jobs to the cluster
with futures.ProcessPoolExecutor(args.num_workers) as workers_pool:
    for idx in range(0,patient_names-load,load):
        jobs.append(
            workers_pool.submit(
                personalization_task, args.file, patient_names[idx:idx+load],args.temp_dir, args.ref, args.output_file
            )
        )
    # push the last file 
    if idx !=len(patient_names):
        workers_pool.submit(
                personalization_task, args.file, patient_names[idx:],args.temp_dir, args.ref, args.output_file
            )
# 5.3 collect the results 
# 5.3.1 initialize counters for success abd failed jobs
num_success=0
num_failure=0
# 5.3.2 wait for the results to finish 
for res in futures.as_completed(jobs):
    try:
        res.result()
        num_success+=1
        print(f"One Batch of patients has been successfully processed, number of successful batches is: {num_success}, total number of batches is :{num_success+num_failure}")
    except:
        num_failure+=1
        print(f"One Batch of patients has not been successfully processed, number of failed batches is: {num_failure}, total number of batches is :{num_success+num_failure}")
# 6. end of execution 
print(f"The file has been processed, the script was called at: {start_time} and it finished at {time.ctime()}, {num_success+num_failure} batches were processed of which {num_failure} failed")






