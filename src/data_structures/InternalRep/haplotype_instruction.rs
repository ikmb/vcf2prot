// load libraries and caret 
use std::{collections::HashMap, panic, usize};
use crate::data_structures::vcf_ds::AltTranscript;
use super::{engines::Engine, task::Task, transcript_instructions::TranscriptInstruction}; 
use rayon::prelude::*; 
use serde::{Deserialize, Serialize};
use crate::data_structures::InternalRep::gir::GIR; 

/// ## Summary
/// An abstraction for a collection on instruction in the same Haplotype of a proband
/// the class derives the Debug, Clone, Serialize and Deserialize traits  
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct HaplotypeInstruction
{
    instructions:Vec<TranscriptInstruction>
}
impl HaplotypeInstruction
{
    /// ## Summary  
    /// Create a new instance using a vector of transcript instruction, i.e. the building block for the instance is a TranscriptInstruction
    /// ##Example 
    /// ```rust
    /// // let's load the modules 
    /// use ppgg_rust::data_structures::{haplotype_instruction::TranscriptInstruction,
    ///                                 transcript_instruction::TranscriptInstruction}
    /// // define an empty vec of transcript instance 
    /// let vec_ins: Vec<TranscriptInstruction> = Vec::new(); 
    /// let hap_ins= HaplotypeInstruction::new(); 
    /// println!("Haplotype instance is: {:#?}",hap_ins); 
    /// ```
    pub fn new(instructions:Vec<TranscriptInstruction>)->Self
    {
        HaplotypeInstruction{instructions}   
    }
    /// ## Summary 
    /// Generate an instance from a vector of AltTranscript, a reference sequence and an execution engine
    pub fn from_vec_t_ins(alt_trans_vec:Vec<AltTranscript>, engine:Engine, ref_seq:&HashMap<String,String>)->Self
    {
        match engine
        {
            Engine::ST=>
            {
                let vec_transcriot_ins= alt_trans_vec.into_iter()
                .map(|alt_transcript| 
                {
                    let name = alt_transcript.name.clone(); 
                    match TranscriptInstruction::from_alt_transcript(alt_transcript, ref_seq)
                    {
                        Ok(res)=>res,
                        Err(err_msg) => TranscriptInstruction::emtpy_t_instruction()
                    }
                })
                .filter(|elem| *elem.get_transcript_name() != "")
                .collect::<Vec<_>>();
                HaplotypeInstruction::new(vec_transcriot_ins)
            }
            Engine::MT | Engine::GPU=>
            {
                let vec_transcriot_ins= alt_trans_vec.into_par_iter()
                .map(|alt_transcript| 
                    {
                        let name = alt_transcript.name.clone(); 
                        match TranscriptInstruction::from_alt_transcript(alt_transcript, ref_seq)
                        {
                            Ok(res)=>res,
                            Err(err_msg) => TranscriptInstruction::emtpy_t_instruction()
                        }
                    })
                    .filter(|elem| *elem.get_transcript_name() != "")
                    .collect::<Vec<_>>();
                HaplotypeInstruction::new(vec_transcriot_ins)
            },
        }
    }
    /// ## Summary
    /// Generate a G Representation from a ref_seq and an execution engine   
    pub fn get_g_rep(&mut self,ref_seq:&HashMap<String,String>, engine:Engine)->GIR
    {
        // Allocate resources 
        let mut results_array=vec!['.'; self.get_size_results_array()];
        let mut alt_array=Vec::with_capacity(self.get_size_alt_array()); 
        let mut reference_array=Vec::with_capacity(self.get_size_ref_array(ref_seq));
        let mut annotation=HashMap::new(); 
        let mut g_rep=Vec::with_capacity(self.get_expected_number_of_tasks());
        // Compute the GIRL representation for each transcript 
        let vec_g_rep= match engine
        {
            Engine::ST=>self.instructions.iter().map(|ins|ins.get_g_rep(ref_seq)).collect::<Vec<_>>(),
            Engine::MT | Engine::GPU =>self.instructions.par_iter().map(|ins|ins.get_g_rep(ref_seq)).collect::<Vec<_>>(),
        };
        // compute some counter 
        let mut ref_counter=0; let mut alt_counter=0; let mut res_counter=0; 
        let mut len_vec=Vec::with_capacity(1000); 
        // loop-and-reindex 
        for g_rep_e in vec_g_rep
        {
            // consume the resources 
            let res=match g_rep_e
            {
                Ok(res)=>res.consumer_and_get_resources(),
                Err(err_msg)=>{println!("While creating instruction for a haplotype, the following error was encountered,{:#?}, skipping this transcript ...\
                Please check your input VCF file, otherwise feel free to contact the developer at: h.elabd@ikmb.uni-kiel.de or at the project webpage: https://github.com/ikmb/ppg", err_msg);continue;},
            };
            // re-index and push the tasks 
            for task in res.0
            {
                len_vec.push(task.get_length());
                g_rep.push(HaplotypeInstruction::update_task(task, &ref_counter, &alt_counter,&res_counter));
            }
            // move the alternative and resource array 
            let (len_alt, len_ref, len_res)=(res.2.len(),res.3.len(),res.4.len()); 
            alt_array.extend(res.2); 
            reference_array.extend(res.3);
            // update the annotation map
            for (key,mut value) in res.1
            {
                value.0+=res_counter;
                value.1+=res_counter;
                annotation.insert(key, value);
            }
            // update the counters 
            ref_counter+=len_ref; 
            alt_counter+=len_alt; 
            res_counter+=len_res; 
        }
        // return the results 
        GIR::new(g_rep, annotation, alt_array, reference_array, results_array)
    }
    /// ## Summary
    /// Update the task index by shifting, i.e. adjusting the position of the task indices 
    fn update_task(mut task:Task,ref_counter:&usize,alt_counter:&usize,res_counter:&usize)->Task
    {
        task=match task.get_execution_stream()
        {
            0=>
            {
                task.shift_start_pos_stream(ref_counter); 
                task
            },
            1=>
            {
                task.shift_start_pos_stream(alt_counter); 
                task
            },
            _=>panic!("Unsupported Stream code: {:#?} from task: {:#?}",task.get_execution_stream(),task)
        }; 
        task.shift_start_pos_res(res_counter);
        task
    }
    /// ## Summary 
    /// compute the size of the results array 
    fn get_size_results_array(&self)->usize
    {
        self.instructions.iter()
        .map(|trans_ins|trans_ins.compute_expected_results_array_size())
        .collect::<Vec<_>>()
        .iter()
        .sum::<usize>()
    }
    /// ## Summary 
    /// compute the size of the alternative array 
    fn get_size_alt_array(&self)->usize
    {
        self.instructions.iter()
        .map(|trans_ins|trans_ins.compute_alt_stream_size())
        .collect::<Vec<_>>()
        .iter()
        .sum::<usize>()
    }
    /// ## Summary 
    /// compute the size of the reference array 
    fn get_size_ref_array(&self, ref_seq:&HashMap<String,String>)->usize
    {
        self.instructions.iter()
        .map(|trans_ins|ref_seq.get(trans_ins.get_transcript_name()).unwrap().len())
        .collect::<Vec<_>>()
        .iter()
        .sum::<usize>()
    }
    /// ## Summary 
    /// compute the number of expected tasks in the instances by summing over the number of tasks in each transcript instruction 
    fn get_expected_number_of_tasks(&self)->usize
    {
        self.instructions.iter()
        .map(|trans_ins|trans_ins.get_num_instructions()*3)
        .collect::<Vec<_>>()
        .iter()
        .sum::<usize>()
    }
}