#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: automatic preperation of samples using a restricted sample size and transcript id
"""
## Load the modules
#------------------
import os 
import pandas as pd
import subprocess as sp
## Define Path to load the data
#------------------------------

OUTPUT_PATH='/work_ifs/sukmb418/human_exome2/predictions_final/fasta'
TEMP_DIR='/work_ifs/sukmb418/human_exome2/predictions_final/fasta'
INPUT_PATH='/work_ifs/sukmb361/2020~2021_HLA_IBD/HLA_IBD_Exome/csq'
LIST_SAMPLES='/work_ifs/sukmb418/human_exome2/sample_sheet.txt'
LIST_TRANS='/work_ifs/sukmb418/human_exome2/predictions_final/ENST_filtered.csv'
REF_PATH='/work_ifs/sukmb418/human_exome2/predictions_final/target_fasta.fasta'

## Load the dataset
#-----------------
txp_list=pd.read_csv(LIST_TRANS,sep=' ').iloc[:,1].to_list()

## GET The Input VCF 
#------------------------------
input_files=[input_file for input_file in os.listdir(INPUT_PATH) if ('bcf.gz' in input_file and 'csi' not in input_file)]

## Create Results
#-------------------------------
for input_file in input_files:
    print(f'Generatating personalized proteomes from: {input_file} ...')
    # create a file to hold the results 
    output_path=os.path.join(OUTPUT_PATH,input_file.split('.')[0])
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
        sp.run(f'module load bcftools; bcftools view -S {LIST_SAMPLES} --threads 16 {os.path.join(INPUT_PATH,input_file)} -O v -o {os.path.join(OUTPUT_PATH,vcf_file_name)}',
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the vcf file failed with the following error: {str(exp)}')
        continue
    # extract the target transcripts  
    try:
        sp.run(f"grep -F '#' {os.path.join(OUTPUT_PATH,vcf_file_name)} > {os.path.join(OUTPUT_PATH,vcf_file_name.strip('.vcf')+'_mini.vcf')}",shell=True,check=True) # grep the header
        sp.run(f"grep -E '{'|'.join(txp_list)}' {os.path.join(OUTPUT_PATH,vcf_file_name)} >> {os.path.join(OUTPUT_PATH,vcf_file_name.strip('.vcf')+'_mini.vcf')}",shell=True,check=True) # extracts the transcripts
    except sp.SubprocessError as exp:
        print(f'Trying to focus on the target transcript from the file cause the following error: {str(exp)}')
    # remove from the major files
    try:
        sp.run(f"rm -rf {os.path.join(OUTPUT_PATH,vcf_file_name)}",shell=True,check=True)
    except sp.SubprocessError as exp:
        print(f'Trying to remove the original file failed with the following error: {str(exp)}')
    # call ppg
    try: 
        sp.run(f"export DEBUG_CPU_EXEC=TRUE; export INSPECT_TXP=TRUE; ./Vcf2prot -f {os.path.join(OUTPUT_PATH,vcf_file_name.strip('.vcf')+'_mini.vcf')} -r {REF_PATH} -o {output_path+'_protoems'} -g mt -vsa &> {os.path.join(TEMP_DIR,vcf_file_name+'Debug_logs')} ", 
            check=True, shell=True)
    except sp.SubprocessError as exp: 
        print(f'Trying to generate the personalized proteomes failed with the following error:{str(exp)}')
        pass 
print(f'Execution finished')
