#!/usr/bin/env python  
"""
@author: Hesham ElAbd
@brief: A simple performance script for benchmarking the vcf2prot against PercisionProDB
@version [script]: 0.2 
@version [vcf2prot]: 0.1.3 
@version [PrecisionProDB]: 1.1.2 
@date [version 0.1] : 11.08.2021
@date [version 0.2]: 20.12.2021
@platform: mac OS/Linux 
"""
## Import the modules
#--------------------  
import os 
import time 
import subprocess as sp 
import pandas as pd 
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    os.system("pip install tqdm")
## define constant parameters 
#---------------------------------------
BASE_DIR='/work_ifs/sukmb418/vcf2prot/'
INPUT_VCF='/work_ifs/sukmb418/vcf2prot/Test_case_chromsome1.vcf'
INPUT_REF='/work_ifs/sukmb418/vcf2prot/References_sequences.fasta'
WRITE_RESULTS_PERCISION_PRO_DB='/work_ifs/sukmb418/vcf2prot/PERCISION_PRO_DB_RESULTS'
VCF2PROT_RES='/work_ifs/sukmb418/vcf2prot/VCF2PROT_RES'
TEMP_WORK_DIR='/work_ifs/sukmb418/vcf2prot/temp_work'
RESULT_PATH='/work_ifs/sukmb418/vcf2prot/benchmark_results'
#---------------------------------------
# PercisionProDB is slow and hence we limit the performance to 128 patient only 
num_patients=[1,2,4,8,16,32,64,128] 
#------------------------------------------
# make output directories to store the results 
#---------------------------------------------
try:
    os.mkdir(WRITE_RESULTS_PERCISION_PRO_DB)
except:pass 

try: 
    os.mkdir(VCF2PROT_RES)
except:pass 

try: 
    os.mkdir(TEMP_WORK_DIR)
except:pass 

try: 
    os.mkdir(RESULT_PATH)
except:pass
## create a logging file
#-----------------------
print(f"Starting the benchmarking loop: time is --> {time.ctime()}")
#------------------------------------------
## allocate some arrays to hold the results 
#------------------------------------------
tool_name,input_size,runtime=[],[],[]
counter=0
for num_pat in tqdm(num_patients):
    # cut the number of patient from the vcf file 
    os.system(f"cut -f 1-9,10-{9+num_pat} {INPUT_VCF} > {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}") # files will be overwritten so we can clean once at the end of the script 
    ## we need to warm up our hard-disk so we run cat on the file and we direct the output into dev null 
    os.system(f" cat {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')} > /dev/null")
    # start the benchmarking by calling vcf2prot 
    ## benchmark vcf2prot        
    start_time=time.time()
    sp.run(f"{os.path.join(BASE_DIR,'vcf2prot')} -f {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}\
            -r {INPUT_REF} -o {VCF2PROT_RES} -g mt -v",
            stdout=sp.DEVNULL,check=True,shell=True) # incase execution failed for whatever reason the whole script shall fail 
    end_time=time.time()
    ## add the result to the list 
    tool_name.append('Vcf2prot')
    input_size.append(num_pat)
    runtime.append(end_time-start_time)
    ## run percision pro database 
    #----------------------------
    samples=['sample'+str(idx+1) for idx in range(num_pat)]
    print(f"Started measuring the performance with precision protein database, starting at: {time.ctime()}")
    start_time=time.time()
    for sample in tqdm(samples):
        sp.run(f"python PrecisionProDB/src/PrecisionProDB.py -m {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}\
             -D GENCODE -o {WRITE_RESULTS_PERCISION_PRO_DB} -s {sample} -g GRCh38.p13.genome.fa.gz -f gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz -p gencode.v39.pc_translations.fa.gz",
        stdout=sp.DEVNULL,check=True,shell=True)
    end_time=time.time()
    ## add the result to the list 
    tool_name.append('PercisionProDB')
    input_size.append(num_pat)
    runtime.append(end_time-start_time)
    counter+=1
    print(f"Current number of iterations is: {counter}, expected number of iterations is: 8, progress is: {(counter/8)*100}%")

# Create a data frame of the results
#-----------------------------------
results_df=pd.DataFrame({
    'tool_name':tool_name,
    'input_size':input_size,
    'runtime':runtime
})
results_df.to_csv(os.path.join(RESULT_PATH,'Performance_Results.tsv'),sep='\t',index=False)
print(f"Execution finished at: {time.ctime()}")
