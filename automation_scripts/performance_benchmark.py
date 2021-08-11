#!/usr/bin/env python  
"""
@author: Hesham ElAbd
@brief: A simple performance script for benchmarking the rust implementation of SIR, ppg
@date: 11.08.2021
"""
## import the modules 
import os 
import time 
import logging
import subprocess as sp 
import pickle
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    os.system("pip install tqdm")
## define constant parameters 
#---------------------------------------
NUM_RUNS_PER_TRIAL=10 
NUM_DISK_WARM_UP=2
INPUT_VCF='/work_ifs/sukmb418/ppg_paper/Test_case_chromsome1.vcf'
INPUT_REF='/work_ifs/sukmb418/ppg_paper/References_sequences.fasta'
TEMP_WORK_DIR='/work_ifs/sukmb418/ppg_paper/temp_work'
TEMP_RESULTS_PATH='/work_ifs/sukmb418/ppg_paper/res'
RESULT_PATH='/work_ifs/sukmb418/ppg_paper/benchmark_results'
#---------------------------------------
## variable parameter
engines=['st','mt','gpu']
num_patients=[1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384]
use_single_thread_write=[True,False]
print(f"Number of runs is: {(NUM_RUNS_PER_TRIAL+NUM_DISK_WARM_UP)*len(engines)*len(num_patients)*len(use_single_thread_write)}")
#------------------------------------------
## create an initialize a dictionary to store the results 
bench_mark_results={  ## loop over the number of patient  
                    elem:{ ## for each trial we loop over engines 
                            eng:{ ## for each engine we loop over wether or not a single thread are used for writting 
                                'single_thread_state'+str(state): []  for state in  use_single_thread_write ## finally we store the number of state in a single thread
                                } for eng in engines 
                            } 
                        for elem in num_patients
                    }
#------------------------------------------
## create a logging file 
logging.basicConfig(filename='runs_log.log',encoding='utf-8',
                    level=logging.DEBUG,format='%(asctime)s %(message)s')
print(f"Starting the benchmarking loop: time is --> {time.ctime()}")
#------------------------------------------
counter=0
for num_pat in tqdm(num_patients):
    # cut the number of patient from the vcf file 
    os.system(f"cut -f 1-9,{9+num_pat} {INPUT_VCF} > {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}") # files will be override so we can clean once at the end
    logging.info(f"Creating the runs input VCF file, the file contains: {num_pat} patients")
    for engine in engines:
        logging.info("Benchmarking with engine {engine}")
        for use_single_write in use_single_thread_write:
            for idx in range(NUM_DISK_WARM_UP+NUM_RUNS_PER_TRIAL): 
                if use_single_write:
                    state_time=time.time()
                    sp.run(f"export DEBUG_GPU=TRUE &&\
                        export INSPECT_TXP=TRUE&&\
                        export INSPECT_INS_GEN=TRUE&&\
                        ./ppg_rust -f {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}\
                             -r {INPUT_REF}, -o {TEMP_RESULTS_PATH} -g {engine} -wv",
                                stdout=sp.DEVNULL,check=True) # incase execution failed for whatever reason the whole script shall fail 
                    end_time=time.time()
                else:
                    state_time=time.time()
                    sp.run(f"export DEBUG_GPU=TRUE &&\
                        export INSPECT_TXP=TRUE&&\
                        export INSPECT_INS_GEN=TRUE&&\
                        ./ppg_rust -f {os.path.join(TEMP_WORK_DIR,f'run_file_with_{num_pat}_patient.vcf')}\
                             -r {INPUT_REF}, -o {TEMP_RESULTS_PATH} -g {engine} -wv",
                                stdout=sp.DEVNULL,check=True)
                    end_time=time.time()
                bench_mark_results[num_pat][engine]['single_thread_state'+str(use_single_write)].append(end_time-state_time)
                counter+=1
                if counter%20==0:
                    print(f"Current number of iterations is: {counter}, expected number of iterations is: 1080, progress is: {(counter/1080)*100}%")
                elif counter%100==0:
                    with open(f"{os.path.join(RESULT_PATH,f'results_at_{counter}_cycle.pickle')}",'wb') as writer_stream:
                        pickle.dump(bench_mark_results,writer_stream)                         
print(f"The benchmark finished at: {time.ctime()}")
print(f"Cleaning temp results and directories ...")
print(f"removing the temp work directory ... starting at {time.ctime()}")
os.system(f"rm -rf {TEMP_WORK_DIR}/*")
print(f"removing the generated fasta files ... starting at  {time.ctime()}")
os.system(f"rm -rf {TEMP_RESULTS_PATH}/*")
print(f"Execution finished at: {time.ctime()}")
