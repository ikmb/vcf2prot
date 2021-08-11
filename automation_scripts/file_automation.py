#!/usr/bin/env python 
"""
@author: Hesham ElAnd
@contact: h.elabd@ikmb.uni-kiel.de 
@brief: A Python script that acts as a wrapper for the underlining executable, i.e. ppg,
This can be used to automate the execution across multiple file where all the input BCF 
files are defined in one directory
@version: 0.0.1 alpha 
@date: 09-08-2021 
"""
## Load the modules
#------------------
import subprocess as sp 
import os 
import argparse
## Define the argument parsers
#-----------------------------
parser=argparse.ArgumentParser(description="A Python script that acts as a wrapper for the underlining executable, i.e. ppg\
 This can be used to automate the execution across multiple file where all the input BCF")

parser.add_argument('--input_path',help="The path to a collection of BCF file that will be sequentially used to generated personalized proteomes.",
                    default='?')
parser.add_argument('--output_path', help="The path to the base directory where results will be written.\
     for each BCF file a directory will be created and the fasta sequence of the proband of the input directory will be written to it.",
    default='?' )
parser.add_argument('--ref_path', help="The path to the reference Fasta file used for generating the results", 
     default='?' )
parser.add_argument('--no_module',help="A boolean flag for whether or not to load the bcftools module, if set the command, module load bcftools will not be executed",
        type='bool',action='store_true')
args=parser.parse_args()
## Validate user input   
#---------------------
input_path_base=args.input_path
if not os.path.exists(input_path_base):
    raise ValueError(f"The input path to files: {input_path_base} is not defined")
output_path_base=args.output_path
if not os.path.exists(input_path_base):
    try: 
        os.mkdir(input_path_base)
    except Exception as exp:
        raise ValueError(f"The provided path to write the write does not exist!!, also creating an exception at the specified path failed with the following error: {exp}")
Reference_name=args.ref_path
if not os.path.exists(input_path_base):
    raise ValueError(f"The provided path for the reference sequence: {Reference_name} does not exists !!")
no_module=args.no_module
## GET NAMES OF THE INPUT FILES
#------------------------------
input_files=[input_file for input_file in os.listdir(input_path_base) if ('bcf.gz' in input_file and 'csi' not in input_file)]
## LOOP OVER THE INPUT FILE AND CALL ppg
#----------------------------------------
for input_file in input_files:
    print(f'Generatating personalized proteomes from: {input_file} ...')
    # create a file to hold the results 
    output_path=os.path.join(output_path_base,input_file.split('.')[0])
    try: 
        os.mkdir(output_path+"_protoems") # make the output file 
    except FileExistsError: 
        pass
    except Exception as exp: 
        print(f'While creating a file to hold the results of file: {input_file} the following error was encounter {str(exp)}')
    # generating the vcf_file
    vcf_file_name=input_file.split('.')[0]+'.vcf'
    # generate the file 
    try:
        if no_module:
            sp.run(f'bcftools view {os.path.join(input_path_base,input_file)} -O v -o {os.path.join(output_path_base,vcf_file_name)}',
            check=True, shell=True)
        else:
            sp.run(f'module load bcftools; bcftools view {os.path.join(input_path_base,input_file)} -O v -o {os.path.join(output_path_base,vcf_file_name)}',
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the vcf file failed with the following error: {str(exp)}')
        continue
    # call ppg
    try: 
        sp.run(f"./ppg_rust -f {os.path.join(output_path_base,vcf_file_name)} -r {Reference_name} -o {output_path+'_protoems'} -g mt -vs &> {os.path.join(output_path_base,vcf_file_name+'Debug_logs')} ", 
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the personalized proteomes failed with the following error:{str(exp)}')
        pass 
    # remove the file: 
    os.system(f'rm -rf {os.path.join(output_path_base,vcf_file_name)}')
print(f'Execution finished')
