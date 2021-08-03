#!/usr/bin/env python 
"""
@author: Hesham ElAnd 
@brief: A Python script that acts as a wrapper for the underlining executable, i.e. PPGG,
This can be used to automate the execution across multiple file 
@version: 0.0.1 alpha 
@date: 02-03-2021 
@copyright: Institute of Clinical Molecular Biology, University of Kiel, Kiel, Germany. 
"""
## LOAD THE MODULES
#------------------
import subprocess as sp 
import os 
## DEFINE THE PATH
#-----------------
input_path_base='/work_ifs/sukmb361/2020~2021_HLA_IBD/HLA_IBD_Exome/csq'
output_path_base='/work_ifs/sukmb418/human_exomes'
Reference_name='References_sequences.fasta'
## GET NAMES OF THE INPUT FILES
#------------------------------
input_files=[input_file for input_file in os.listdir(input_path_base) if ('bcf.gz' in input_file and 'csi' not in input_file)]
## LOOP OVER THE INPUT FILE AND CALL PPGG
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
        sp.run(f'module load bcftools; bcftools view {os.path.join(input_path_base,input_file)} -O v -o {os.path.join(output_path_base,vcf_file_name)}',
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the vcf file failed with the following error: {str(exp)}')
        continue
    # call PPGG
    try: 
        sp.run(f"./ppgg_rust -f {os.path.join(output_path_base,vcf_file_name)} -r {Reference_name} -o {output_path+'_protoems'} -g mt -vs &> {os.path.join(output_path_base,vcf_file_name+'Debug_logs')} ", 
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the personalized proteomes failed with the following error:{str(exp)}')
        pass 
    # remove the file: 
    os.system(f'rm -rf {os.path.join(output_path_base,vcf_file_name)}')
print(f'Execution finished')
