use std::collections::HashMap;
use std::collections::HashSet;
use std::panic;
use std::usize;
use crate::data_structures::InternalRep::gir; 
use crate::data_structures::InternalRep::instruction;
use crate::data_structures::vcf_ds; 
use crate::data_structures::InternalRep::task::Task;
use serde::{Deserialize, Serialize};
use std::iter::FromIterator;
use super::instruction::Instruction;

/// A representation for a collection of mutation in a transcript, where mutations have been already encoded into instructions 
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct TranscriptInstruction
{
    transcript_name:String,
    ref_len:usize, 
    instructions:Vec<instruction::Instruction>
}
impl TranscriptInstruction
{
    /// ## Summary 
    /// Create a new instruction from a transcript name and vector of instruction containing 
    /// the mutation in the transcript and the reference length
    pub fn new(transcript_name:String, ref_len:usize, instructions:Vec<instruction::Instruction>)->Self
    {
        TranscriptInstruction{transcript_name,ref_len,instructions}
    }
    /// ## Summary 
    /// Create a new instruction from a transcript name and vector of instruction containing 
    /// the mutation in the transcript and the reference length
    pub fn empty_t_instruction()->Self
    {
        let ins_vec:Vec<instruction::Instruction>=Vec::new(); 
        TranscriptInstruction{transcript_name:"".to_string(),ref_len:1,instructions:ins_vec}
    }
    /// ## Summary 
    /// Create a new instance from the alt-Transcript instance along with a reference hashmap 
    /// of sequence names 
    pub fn from_alt_transcript(mut alt_transcript:vcf_ds::AltTranscript, ref_seqs:&HashMap<String,String>)->Result<Self,String>
    {
        alt_transcript.sort_alterations();// sort alteration 
        let transcript_name=alt_transcript.name.clone(); 
        let ref_len=match ref_seqs.get(&transcript_name)
        {
            Some(sequence)=>sequence.len(),
            None=>return Err(format!("The provided transcript name: {} is not in the reference sequence", &transcript_name))
        };
        let mut instructions= Vec::with_capacity(alt_transcript.alts.len()); 
        let mut ins_idx=Vec::with_capacity(alt_transcript.alts.len());
        for mutation in alt_transcript.alts.iter()
        {
            let instruction=instruction::Instruction::from_mutation(mutation,&alt_transcript.alts);
            if instruction.get_code()!='E'
            {
                ins_idx.push(instruction.get_position_ref());
                instructions.push(instruction)
            }
        }
        if instructions.len()==0
        {
            return Err(format!("The provided transcript name: {} has {} mutations none of them is supported, skipping this transcript", &transcript_name,alt_transcript.alts.len()))
        }
        match std::env::var("INSPECT_INS_GEN")
        {
            Ok(_)=>
            {
                println!("Inspect the generation of instructions from transcript: {}",&transcript_name); 
                println!("Checking for indexing errors ...."); 
                let index = HashSet::<usize>::from_iter(ins_idx.into_iter());
                // check that all instructions have unique start points 
                //-----------------------------------------------------
                if index.len()!=instructions.len()
                {
                    match std::env::var("PANIC_INSPECT_ERR")
                    {
                        Ok(_)=>
                        {
                            panic!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                            Nevertheless, some mutations at the same position were observed, skipping this mutations ...
                            ",transcript_name,&alt_transcript.alts,&instructions);   
                        },
                        Err(_)=>
                        {
                            println!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                            Nevertheless, some mutations at the same position were observed, skipping this mutations ...
                            ",transcript_name,&alt_transcript.alts,&instructions);   
                            return Err("Encountered two mutation at the same time".to_string()); 
                        }  
                    }
                }
                println!("Checking for mutational overlaps ...");
                // generate the pair of adjacent instructions
                //-------------------------------------------
                if instructions.len()>1
                {
                    let consequent_pairs=instructions[..instructions.len()-1]
                            .iter()
                            .enumerate()
                            .map(|(idx,ins)|(ins,&instructions[idx+1]))
                            .collect::<Vec<_>>(); 
                    for (ins1,ins2) in consequent_pairs
                    {
                        if ins2.get_position_res()<=(ins1.get_position_res()+ins1.get_data().len()+1)
                        {
                            match std::env::var("PANIC_INSPECT_ERR")
                            {
                                Ok(_)=>
                                {
                                    panic!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                                    Nevertheless, some mutations overlap, skipping this mutations ...
                                    ",transcript_name,&alt_transcript.alts,&instructions);   
                                },
                                Err(_)=>
                                {
                                    println!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                                    Nevertheless, some mutations overlap, skipping this mutations ...
                                    ",transcript_name,&alt_transcript.alts,&instructions);   
                                    return Err("Encountered two mutation at the same time".to_string()); 
                                }  
                            }
                        }
                        // add a conditional state for handling inframe_deletion 
                        if ins1.get_code()=='C' || ins1.get_code()=='D'
                        {
                            if ins2.get_position_ref()<=(ins1.get_position_res()+ins1.get_length()+1) // for example a deletion at position 84 with length 5 (ins1) and ins2 is a mutation at site 88 
                            {
                                match std::env::var("PANIC_INSPECT_ERR")
                                {
                                    Ok(_)=>
                                    {
                                        panic!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                                        Nevertheless, some mutations overlap, skipping this mutations ...
                                        ",transcript_name,&alt_transcript.alts,&instructions);   
                                    },
                                    Err(_)=>
                                    {
                                        println!("Critical error was encountered: for transcript: {} the following mutations were observed, {:#?} which translate to the following instructions: {:#?}\
                                        Nevertheless, some mutations overlap, skipping this mutations ...
                                        ",transcript_name,&alt_transcript.alts,&instructions);   
                                        return Err("Encountered two mutation at the same time".to_string()); 
                                    }  
                                }
                            }
                        }                        
                    }      
                }
            }
            Err(_)=>()
        }
        match std::env::var("DEBUG_TXP")
        {
            Ok(txp_id)=>
            {
                if transcript_name==txp_id
                {
                    println!("Current mutations are : {:#?} and it has been translated into the following instruction: {:#?}", &alt_transcript,&instructions)
                }
            },
            Err(_)=>()
        }
        Ok(TranscriptInstruction::new(transcript_name,ref_len,instructions))
    }
    /// ## Summary 
    /// Return the number of instruction in the transcript 
    pub fn get_num_instructions(&self)->usize
    {
        self.instructions.len()
    }
    /// ## Summary 
    /// Return the transcript name 
    pub fn get_transcript_name(&self)->&String
    {
        &self.transcript_name
    }
    /// ## Summary 
    /// Return a mutable reference to the instance vector of instructions 
    pub fn get_mut_instruction(&mut self)->&mut Vec<Instruction>
    {
        &mut self.instructions
    }
    /// ## Summary 
    /// Compute the size, i.e. the number of chars, in the alternative stream array 
    pub fn compute_alt_stream_size(&self)->usize
    {
        let mut counter=0; 
        for ins in self.instructions.iter()
        {
            counter+=ins.get_data().len(); 
        }
        counter
    }
    /// Inspect the TranscriptInstruction instance and compute the size of the results array
    /// ## Experiment 
    ///```  
    /// use ppg_rust::data_structures::InternalRep::transcript_instructions::TranscriptInstruction; 
    /// use ppg_rust::data_structures::mutation_ds::Mutation; 
    /// use ppg_rust::data_structures::InternalRep::instruction;
    /// let test_case=vec!["frameshift".to_string(), "ENST00000510017".to_string(), "40VGLHFWTM*>40VDSTFGQC".to_string()];
    /// let test_mutation=Mutation::new(test_case).unwrap();
    /// let ins=instruction::Instruction::from_mutation(&test_mutation); 
    /// let mut ins_vec=Vec::with_capacity(2);
    /// ins_vec.push(ins);
    /// let test_alt_transcript=TranscriptInstruction::new("Test1".to_string(), 50, ins_vec);
    /// println!("{:#?}",test_alt_transcript);
    /// assert_eq!(test_alt_transcript.compute_expected_results_array_size(),47);
    ///```  
    pub fn compute_expected_results_array_size(&self)->usize
    {
        let mut expected_size=0; 
        for ins in self.instructions.iter()
        {
            match ins.get_code()
            {
                'U' | '0'=> {expected_size -= self.ref_len as i32;break;},
                'F' => expected_size += ins.get_data().len()  as i32 - (self.ref_len as i32 - ins.get_position_ref() as i32), 
                'R' => 
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();

                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size += ins.get_data().len()  as i32 - (self.ref_len as i32 - ins.get_position_ref() as i32)
                    }
                },
                'G' | 'X' => expected_size -= self.ref_len as i32 - ins.get_position_ref() as i32, 
                'M' | 'N' | '2' => (), 
                'L' =>
                {
                    if ins.get_position_ref()+1 as usize == self.ref_len ||  ins.get_position_ref() == self.ref_len
                    {
                        expected_size += ins.get_data().len() as i32
                    }
                    else
                    {
                        expected_size += ins.get_data().len()  as i32 - (self.ref_len as i32 - ins.get_position_ref() as i32)
                    }
                },
                'I' => expected_size += ins.get_data().len() as i32 -1 as i32, // e.g. 125Y>125YRR,
                'J' => 
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();
                    // check if there is other mutation before hand 
                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size += ins.get_data().len() as i32 -1 as i32;
                    }
                }
                'D' => expected_size -= ins.get_length() as i32, 
                'C' => 
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();
                    // check if there is other mutation before hand 
                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size -= ins.get_length() as i32;
                    }
                },
                'K' => 
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();
                    // check if there is other mutation before hand 
                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size+= ins.get_data().len()  as i32 - (self.ref_len as i32 -ins.get_position_ref() as i32); 
                    }
                },
                'Q' => 
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();
                    // check if there is other mutation before hand 
                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size+= ins.get_data().len()  as i32 - (self.ref_len as i32 -ins.get_position_ref() as i32); 
                    }
                },
                'A' =>
                {
                    // get the index of the current instruction
                    let index=self.instructions.iter().position(|i|i==ins).unwrap();
                    // check if there is other mutation before hand 
                    if !self.instructions[..index].iter().any(|ins|ins.get_code()=='G' || ins.get_code()=='F')
                    {
                        expected_size-= self.ref_len as i32 - ins.get_position_ref() as i32; 
                    }
                },
                'B' => expected_size-= self.ref_len as i32 - ins.get_position_ref() as i32 - ins.get_length() as i32,
                'P' => expected_size-= ins.get_length() as i32 ,
                'Z' => (),
                'T' => expected_size-= self.ref_len as i32 - ins.get_position_ref() as i32,
                'W' => expected_size+= ins.get_data().len() as i32,
                'Y' => expected_size+= ins.get_data().len()  as i32 - (self.ref_len as i32 -ins.get_position_ref() as i32)  +1 as i32, 
                '3' => expected_size+= ins.get_data().len() as i32 - ins.get_length() as i32 ,
                _=>panic!("instruction: {:#?} is not supported", ins),
            }
        }
        let size = (self.ref_len as i32 + expected_size) as usize;
        match std::env::var("DEBUG_TXP")
        {
            Ok(txp_id)=>
            {
                if self.transcript_name==txp_id
                {
                    println!("Transcript Name is: {}, Input protein length is: {}, mutations are: {:?} and computed results is: {}", self.transcript_name,self.ref_len, self.instructions,size);
                }
            },
            Err(_)=>()
        }
        size
    }
    /// Return an GIR  of the instances 
    /// ## Example
    ///```  
    /// let name="ENST00000406869".to_string(); 
    /// let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|10V>10H|1936821C>T".to_string()]; 
    /// let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
    /// println!("{:#?}",alt_transcript); 
    /// let mut reference=HashMap::new(); 
    /// reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
    /// let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
    /// let test_gir=res.get_g_rep(&reference); 
    /// println!("{:#?}",test_gir); 
    ///```
    pub fn get_g_rep(&self, ref_seqs:&HashMap<String,String>)->Result<gir::GIR,String>
    {        
        // handle the case with start-lost and 'U' code
        if self.instructions.iter().any(|ins| ins.get_code()=='0' || ins.get_code()=='U') || self.instructions.len() ==0
        {
            let mut annotations=HashMap::new();
            annotations.insert(self.transcript_name.clone(), (0 as usize,0 as usize));
            return Ok(gir::GIR::new(Vec::new(),annotations,Vec::new(),Vec::new(),Vec::new())); 
        }
        // allocate arrays:
        //-----------------
        let mut vec_tasks=Vec::with_capacity(2*self.instructions.len()); 
        let mut alt_array=Vec::with_capacity(self.compute_alt_stream_size());
        let res_array=vec!['.'; self.compute_expected_results_array_size()];
        let ref_stream=ref_seqs.get(&self.transcript_name).unwrap().chars().collect::<Vec<char>>();
        // push the instruction 
        //---------------------
        // base instruction
        vec_tasks.push(TranscriptInstruction::build_base_instruction(&self.instructions[0],&self.ref_len)); 
        // loop over all instructions
        for ins in self.instructions.iter()
        {
            let (task1, task2)=match TranscriptInstruction::to_task(ins, &self.instructions,
                &mut alt_array, &vec_tasks, ref_stream.len())
            {
                Ok(res)=>res,
                Err(err)=> return Err(format!("Translating {} failed with the following error: {:?}",self.transcript_name, err))
            };
            if !(*task1.get_execution_stream() == 2 as u8)
            {
                vec_tasks.push(task1);
            }
            if !(*task2.get_execution_stream() == 2 as u8)
            {
                vec_tasks.push(task2);
            }           
        }
        match std::env::var("DEBUG_TXP")
        {
            Ok(txp_id)=>
            {
                if self.transcript_name==txp_id
                {
                    println!("Vector of tasks is: {:#?}",vec_tasks);
                }
            },
            Err(_)=>()
        }
        // add the instruction to the array 
        let mut annotations=HashMap::new();
        annotations.insert(self.transcript_name.clone(), (0  as usize, self.compute_expected_results_array_size())); 
        match std::env::var("INSPECT_TXP")
        {
            Ok(_)=>
            {
                let expected_size= res_array.len();
                let mut counter=0; 
                println!("Validating the generated instructions for transcript: {} ....", self.transcript_name);
                for idx in 1..vec_tasks.len()
                {
                    if vec_tasks[idx].get_start_pos_res()!=vec_tasks[idx-1].get_start_pos_res() + vec_tasks[idx-1].get_length()
                    {
                        println!("************ CPU Execution Table *********");
                        println!("index\tstream\tstart_position\tlength\tposition_results\t");
                        for idx in 0..vec_tasks.len()
                        {
                            println!(
                                "{}\t{}\t{}\t{}\t{}\t",idx,vec_tasks[idx].get_execution_stream(),
                                vec_tasks[idx].get_start_pos(),
                                vec_tasks[idx].get_length(),
                                vec_tasks[idx].get_start_pos_res()
                                );
                        }
                        panic!("Critical failure in the calculations was encountered: position: {} the sum {} does not equal previous inputs: {} and {} \n",
                        idx,vec_tasks[idx].get_start_pos_res(),vec_tasks[idx-1].get_start_pos_res(),vec_tasks[idx-1].get_length());
                    }
                    counter+=vec_tasks[idx].get_length()
                }
                counter+=vec_tasks[0].get_length(); 
                if expected_size != counter
                {
                    panic!("In transcript: {}, size mismatched between the expected size: {} and actual size: {}.",
                    self.transcript_name, expected_size,counter); 
                }
            },
            Err(_)=>()
        }

        Ok(gir::GIR::new(vec_tasks, annotations,alt_array,ref_stream,res_array))
    }
    /// Takes an instruction and returns two tasks, the first is the  execution task for the instruction
    /// and the second is the taskdescribe the copying of the reference untill the end of the sequence or
    /// untill the next instruction.
    /// ## Example
    ///```
    /// let name="ENST00000406869".to_string(); 
    /// let mutations=vec![
    ///    "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|5G>5GTEST|1936821C>T".to_string(),
    ///    "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|10V>10VECT|1936821C>T".to_string(),
    /// ]; 
    /// let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
    /// println!("{:#?}",alt_transcript); 
    /// let mut reference=HashMap::new(); 
    /// reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
    /// let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
    /// let test_gir=res.get_g_rep(&reference); 
    /// println!("{:#?}",test_gir); 
    /// let (res_array, _)=test_gir.execute(Engine::ST);
    /// let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
    /// println!("Input Sequence is:  ==>{:#?}",&ref_string);
    /// let res_string=res_array.iter().collect::<String>();
    /// println!("Result sequence is: ==>{:#?}",res_string);
    /// assert_eq!(ref_string.len()+ 7 as usize, res_string.len());
    ///```
    fn to_task(instruction:&instruction::Instruction, vec_instruction:&Vec<instruction::Instruction>, 
                    alt_stream: &mut Vec<char>, vec_tasks: &Vec<Task>, ref_len:usize)->Result<(Task,Task),String>
    {
        
        
        let ins_task= match instruction.get_code() 
        {
            'M' => TranscriptInstruction::get_task_from_missense(instruction, alt_stream, vec_tasks),
            'N' => TranscriptInstruction::get_task_from_missense(instruction, alt_stream, vec_tasks),
            'F' => TranscriptInstruction::get_task_from_frameshift(instruction, alt_stream, vec_tasks),
            'R' => TranscriptInstruction::get_task_from_frameshift(instruction, alt_stream, vec_tasks),
            'G' => TranscriptInstruction::get_task_from_stop_gained(instruction, alt_stream, vec_tasks),
            'X' => TranscriptInstruction::get_task_from_s_stop_gained(instruction,alt_stream, vec_tasks),
            'L' => TranscriptInstruction::get_task_from_stop_lost(instruction, alt_stream, vec_tasks),
            'I' => TranscriptInstruction::get_task_from_inframe_insertion(instruction, alt_stream, vec_tasks), 
            'J' => TranscriptInstruction::get_task_from_inframe_insertion(instruction, alt_stream, vec_tasks), 
            'D' => TranscriptInstruction::get_task_from_inframe_deletion(instruction, alt_stream, vec_tasks),
            'C' => TranscriptInstruction::get_task_from_inframe_deletion(instruction, alt_stream, vec_tasks),
            'K' => TranscriptInstruction::get_task_from_frameshift(instruction,alt_stream, vec_tasks),
            'Q' | 'Z' | 'P'=> Task::new(2, 0, 0,0), // a phi-instruction ,
            'A' => TranscriptInstruction::get_task_from_stop_gained(instruction,alt_stream, vec_tasks),
            'B' => TranscriptInstruction::get_task_from_frameshift(instruction,alt_stream, vec_tasks),
            'T' => TranscriptInstruction::get_task_from_stop_gained(instruction,alt_stream, vec_tasks),
            'W' => TranscriptInstruction::get_task_from_stop_lost(instruction,alt_stream, vec_tasks),
            'Y' => TranscriptInstruction::get_task_from_frameshift(instruction,alt_stream, vec_tasks),
            '2' => TranscriptInstruction::get_task_from_instruction_2(instruction, alt_stream, vec_tasks),
            '3' => TranscriptInstruction::get_task_from_instruction_3(instruction, alt_stream, vec_tasks),
            _=>panic!("Instruction: {:#?} is not supported",instruction)
        };
        let last_ins=vec_instruction.last().unwrap() == instruction; 
        let last_ins= match  last_ins
        {
            true =>
            {
                let last_task_type=['K','Y','Q','A','B','P','Z','T','W','Z','T','W','Z','G','F','R','L','X'].iter().any(|c|*c==instruction.get_code()); 
                
                match last_task_type
                {
                    true => Task::new(2, 0,0,0),
                    false => TranscriptInstruction::add_last_instruction(ref_len, instruction,ins_task.get_start_pos_res()+ins_task.get_length())
                }
            },
            false =>
            {
                let last_task_type=['K','Q','A','B','P','Z','T','W','Z','T','W','Z','G','F','R','L'].iter().any(|c|*c==instruction.get_code()); 
                match last_task_type
                {
                    true => return Err(format!("Instruction: {:?} must be the last mutation in a transcript, current vector of instructions is: {:?}",instruction,vec_instruction)),
                    false => TranscriptInstruction::add_till_next_ins(instruction, vec_instruction,&ins_task, &ref_len)
                }
            }
        };        
        Ok((ins_task,last_ins))   
    }
    /// ## Summary 
    /// add an instruction to copy the data until the next instruction in the vector of mutations
    fn add_till_next_ins(ins:&instruction::Instruction, instructions:&Vec<instruction::Instruction>, last_task:&Task, ref_len:&usize)->Task
    {
        let position=instructions.iter().position(|inst_cmp|inst_cmp==ins).unwrap();
        let next_ins=&instructions[position+1 as usize];
        match ins.get_code()
        {
            'D' | 'C'=>
            {
                if next_ins.get_position_ref() == ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else if (ins.get_position_ref() + ins.get_length())==next_ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else if next_ins.get_code()=='L' && next_ins.get_position_ref()+1 as usize == *ref_len && (ins.get_position_ref()+ins.get_length() +1 as usize)==next_ins.get_position_ref() 
                {
                    let stat_pos=ins.get_position_ref()+ins.get_length() +1 as usize; 
                    return Task::new(0, stat_pos, 1,last_task.get_start_pos_res()+last_task.get_length())
            
                }
                else 
                {
                    let stat_pos=ins.get_position_ref()+ins.get_length() +1 as usize; 
                    let res= panic::catch_unwind(||
                    {
                        Task::new(0, stat_pos, 
                            next_ins.get_position_ref()- stat_pos, 
                            last_task.get_start_pos_res()+last_task.get_length())
                    }
                    );
                    match res
                    {
                        Ok(res)=>res, 
                        Err(err)=>panic!("This error was encountered: {:#?} while generating a task from: {:#?}, {:#?}, {:#?},{:#?},\
                        The current instruction is: {:#?} while all instructions are: {:#?}",
                            err, stat_pos,next_ins.get_position_ref(),last_task.get_start_pos_res(),last_task.get_length(), ins,instructions )
                    }   
                }
            }, 
            '2'  =>
            {
                if next_ins.get_position_ref() == ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else if (ins.get_position_ref() + ins.get_length())==next_ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else
                {
                    let start_pos=ins.get_position_ref() + ins.get_length();
                    Task::new(0, 
                        start_pos, 
                    next_ins.get_position_ref() - start_pos, 
                    last_task.get_start_pos_res()+last_task.get_length())
                }
            },
            '3' =>
            {
                if next_ins.get_position_ref() == ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else if (ins.get_position_ref() + ins.get_length())==next_ins.get_position_ref()
                {
                    return Task::new(2, 0, 0,0)
                }
                else
                {
                    let start_pos=ins.get_position_ref() + ins.get_length(); 
                    Task::new(0, 
                        start_pos, 
                    next_ins.get_position_ref() - start_pos, 
                    last_task.get_start_pos_res()+last_task.get_length())
                }
            },
            _=>
            {
                if next_ins.get_position_ref() == ins.get_position_ref()
                {
                    Task::new(2, 0, 0,0)
                }
                else
                {
                    if next_ins.get_code()=='L'
                    {
                        if next_ins.get_position_ref()+1 as usize ==*ref_len
                        {
                            return Task::new(0, ins.get_position_ref()+1 as usize,
                        next_ins.get_position_ref() - ins.get_position_ref(),
                        last_task.get_start_pos_res()+last_task.get_length())
                        }
                        else
                        {
                            return Task::new(0, ins.get_position_ref()+1 as usize,
                        next_ins.get_position_ref() - (1 as usize) - ins.get_position_ref(),
                        last_task.get_start_pos_res()+last_task.get_length())
                        }
                    }
                    let res=panic::catch_unwind( ||
                    {
                        Task::new(0, ins.get_position_ref()+1 as usize,
                    next_ins.get_position_ref() - (1 as usize) - ins.get_position_ref(),
                    last_task.get_start_pos_res()+last_task.get_length())
                    }); 
                    match res {
                        Ok(task)=>task,
                        Err(err_msg)=>
                        {
                            println!("Instructions are : {:#?}",instructions); 
                            panic!("The following error: {:#?} cause by: {}, {}",
                        err_msg,next_ins.get_position_ref(), ins.get_position_ref())
                        }
                    }
                    
                }
            }
        }
    }
    /// ## Summary 
    /// Add the last instruction, i.e. adds an instruction that describing copying from the last mutation instruction 
    /// in the transcript until the end of the transcript 
    fn add_last_instruction(ref_seq:usize, instruction:&instruction::Instruction, pos_res_array:usize)->Task
    {
        match instruction.get_code()
        {
            'D' | 'C'=>
            {
                Task::new(0, instruction.get_position_ref()+instruction.get_length()+ 1 as usize, 
                ref_seq-instruction.get_position_ref()-instruction.get_length() -1 as usize, 
                pos_res_array)
            },
            '2' | '3' =>
            {
                Task::new(0, instruction.get_position_ref()+instruction.get_length(), 
            ref_seq-instruction.get_position_ref()-instruction.get_length() as usize, pos_res_array)
            },
            _=> Task::new(0, instruction.get_position_ref()+1, 
            ref_seq-instruction.get_position_ref()-1 as usize, pos_res_array)
        }        
    }
    /// ## Summary 
    /// returns a Task from a missense mutation encoded as an instruction 
    fn get_task_from_missense(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length() ;
        alt_stream.extend(instruction.get_data().iter());
        let pos_altstream= match alt_stream.len()
        {
            0=>0,
            _=>alt_stream.len()-1
        }; 
        Task::new(1, pos_altstream, 1,pos_result)
    }
    /// ## Summary 
    /// returns a Task from a frameshift mutation encoded as an instruction 
    fn get_task_from_frameshift(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= match alt_stream.len()
        {
            0=>0,
            _=>alt_stream.len()-1
        }; 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length() ;
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream, instruction.get_length(),pos_result)

    }
    /// ## Summary 
    /// returns a Task from a stop-gained mutation encoded as an instruction 
    fn get_task_from_stop_gained(_instruction:&instruction::Instruction, _alt_stream:&mut Vec<char>,
        _vec_tasks:&Vec<Task>)->Task
    {
        Task::new(2, 0, 0,0)
    }
    /// ## Summary 
    /// returns a Task from an *stop-gained mutation encoded as an instruction 
    fn get_task_from_s_stop_gained(_instruction:&instruction::Instruction, _alt_stream:&mut Vec<char>,
        _vec_tasks:&Vec<Task>)->Task
    {
        Task::new(2, 0, 0,0)
    }
    /// ## Summary 
    /// returns a Task from a stop_lost mutation encoded as an instruction 
    fn get_task_from_stop_lost(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= match alt_stream.len()
        {
            0=>0,
            _=>alt_stream.len()-1
        }; 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res()+ last_task.get_length() ;
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream, instruction.get_data().len(),pos_result)
    }
    /// ## Summary 
    /// Build the base instruction for a transcript where the reference is copied until the first instruction 
    fn build_base_instruction(instruction:&instruction::Instruction, ref_len:&usize)->Task
    {
        match instruction.get_code()
        {
            'Z' | 'Y'  => Task::new(0, 0, instruction.get_position_ref()+1, 0),
            'L'=>
            {
                if instruction.get_position_ref()+1 as usize ==*ref_len
                {
                    Task::new(0, 0, instruction.get_position_ref()+1, 0)
                }
                // edit here 
                else if instruction.get_position_ref()==*ref_len
                {
                    Task::new(0, 0, instruction.get_position_ref(), 0)
                }
                else 
                {
                    Task::new(0, 0, instruction.get_position_res(), 0)
                }
            },
            _=> Task::new(0, 0, instruction.get_position_ref(), 0)
        }
    }
    /// ## Summary 
    /// Returns a Task from an inframe_insertion mutation encoded as an instruction 
    fn get_task_from_inframe_insertion(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= alt_stream.len(); 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length() ;
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream, instruction.get_length(),pos_result)
    }
    /// ## Summary 
    /// Returns a Task from an inframe_deletion mutation encoded as an instruction 
    fn get_task_from_inframe_deletion(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= alt_stream.len(); 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length();
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream,  instruction.get_data().len(),pos_result)
    }
    /// ## Summary 
    /// Returns a Task from a 2 instruction, see the instruction module for more details   
    fn get_task_from_instruction_2(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= alt_stream.len(); 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length();
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream, instruction.get_length(),pos_result)
    }
    /// ## Summary 
    /// Returns a Task from a 3 instruction, see the instruction module for more details 
    fn get_task_from_instruction_3(instruction:&instruction::Instruction, alt_stream:&mut Vec<char>,
        vec_tasks:&Vec<Task>)->Task
    {
        let pos_altstream= alt_stream.len(); 
        let last_task=vec_tasks.last().unwrap(); 
        let pos_result=last_task.get_start_pos_res() + last_task.get_length();
        alt_stream.extend(instruction.get_data().iter());
        Task::new(1, pos_altstream, instruction.get_data().len(),pos_result)
    }
}

#[cfg(test)]
pub mod test_transcript_instruction
{
    use super::*; 
    use crate::data_structures::mutation_ds::Mutation; 
    use super::super::engines::Engine; 
    #[test]
    pub fn test_expected_result_array_length()
    {
        let test_case=vec!["frameshift".to_string(),"ENST00000510017".to_string(), "40VGLHFWTM*>40VDSTFGQC".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=instruction::Instruction::from_mutation(&test_mutation); 
        let mut ins_vec=Vec::with_capacity(2);
        ins_vec.push(ins);
        let test_alt_transcript=TranscriptInstruction::new("Test1".to_string(), 50, ins_vec);
        println!("{:#?}",test_alt_transcript);
        assert_eq!(test_alt_transcript.compute_expected_results_array_size(),47);
    }
    #[test]
    pub fn test_get_task_from_frameshift()
    {
        let test_case=vec!["frameshift".to_string(),"ENST00000510017".to_string(), "40VGLHFWTM*>40VDSTFGQC".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=instruction::Instruction::from_mutation(&test_mutation); 
        let mut vec_tasks=Vec::with_capacity(2);
        vec_tasks.push(Task::new(0, 1, 15, 15));
        let mut alt_stream=Vec::with_capacity(100);
        let task=TranscriptInstruction::get_task_from_frameshift(&ins,&mut alt_stream, &mut vec_tasks);
        println!("Defining Task: {:#?}",Task::new(0, 1, 15, 15));
        println!("Defining second task: {:#?}",&task);
        assert_eq!(task,Task::new(1, 0, 8, 30));
    }
    #[test]
    pub fn test_get_task_from_stop_gained()
    {
        let test_case=vec!["stop_gained".to_string(),"ENST00000510017".to_string(), "40VGLHFWTM*>40*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=instruction::Instruction::from_mutation(&test_mutation); 
        let mut vec_taks=Vec::with_capacity(2);
        vec_taks.push(Task::new(0, 0, 39, 0));
        let mut alt_stream=Vec::with_capacity(100);
        let task=TranscriptInstruction::get_task_from_stop_gained(&ins,&mut alt_stream,&mut vec_taks);
        println!("Base Task: {:#?}",Task::new(0, 0, 40, 0));
        println!("Defining second task: {:#?}",&task);
        assert_eq!(task,Task::new(2, 0, 0, 0));
    }
    #[test]
    pub fn test_get_task_from_stop_lost() 
    {
        let test_case=vec!["stop_lost".to_string(),"ENST00000650310".to_string(), "489*>489S".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=instruction::Instruction::from_mutation(&test_mutation); 
        let mut vec_taks=Vec::with_capacity(2);
        vec_taks.push(Task::new(0, 0, 488, 0));
        let mut alt_stream=Vec::with_capacity(100);
        let task=TranscriptInstruction::get_task_from_stop_lost(&ins, &mut alt_stream,&mut vec_taks);
        println!("Base Task: {:#?}",Task::new(0, 0, 488, 0));
        println!("Defining second task: {:#?}",&task);
        assert_eq!(task,Task::new(1, 0, 1, 488));
        assert_eq!(alt_stream.len(),1);
        assert_eq!(alt_stream[0],'S');   
    }
    #[test]
    pub fn test_get_task_from_inframe_insersion() 
    {
        let test_case=vec!["inframe_insertion".to_string(),"ENST00000484547".to_string(), "125Y>125YRR".to_string()];        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=instruction::Instruction::from_mutation(&test_mutation); 
        let mut vec_taks=Vec::with_capacity(2);
        vec_taks.push(Task::new(0, 0, 124, 0));
        let mut alt_stream=Vec::with_capacity(10);
        let task=TranscriptInstruction::get_task_from_stop_lost(&ins, &mut alt_stream,&mut vec_taks);
        println!("Base Task: {:#?}",Task::new(0, 0, 124, 0));
        println!("Defining second task: {:#?}",&task);
        assert_eq!(task,Task::new(1, 0, 3, 124));
        assert_eq!(alt_stream.len(),3);
        assert_eq!(alt_stream[0],'Y');   
        assert_eq!(alt_stream[1],'R');   
        assert_eq!(alt_stream[2],'R');   
    }
    #[test]
    pub fn test_correct_translation_1()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|5G>5H|1936821C>T".to_string()]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, res_map)=test_gir.execute(Engine::ST);
        println!("Res");
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        test_equal_expect(&ref_string,&res_string,vec![4;1]); 
        println!("Result sequence is: ==>{:#?}",res_string);
        
    }
    fn test_equal_expect(string1:&String, string2:&String, indices:Vec<usize>)
    {
        assert_eq!(string1.len(),string2.len());
        let char_vec1=string1.chars().collect::<Vec<char>>();
        let char_vec2=string2.chars().collect::<Vec<char>>();
        for idx in 0..string1.len()
        {
            if !(indices.iter().any(|num|idx==*num))
            {
                assert_eq!(char_vec1[idx],char_vec2[idx]);
            }
        }
    }
    #[test]
    pub fn test_correct_translation_2()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|10V>10H|1936821C>T".to_string()]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, res_map)=test_gir.execute(Engine::ST);
        println!("Res");
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        test_equal_expect(&ref_string,&res_string,vec![9;1]); 
    }
    #[test]
    pub fn test_correct_translation_3()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|10V>10H|1936821C>T".to_string(),
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|20F>20K|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(&ref_string.len(),&res_string.len());
        test_equal_expect(&ref_string,&res_string,vec![9,19]); 
    }
    #[test]
    fn test_correct_translation_4()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|5G>5T|1936821C>T".to_string(),
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|10V>10E|1936821C>T".to_string(),
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|15R>15S|1936821C>T".to_string(),
            "*missense|MAD1L1|ENST00000406869|protein_coding|-|20F>20T|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(&ref_string.len(),&res_string.len());
        test_equal_expect(&ref_string,&res_string,vec![4,9,14,19]); 
    }
    #[test]
    fn test_correct_translation_5()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|5G>5GTEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(ref_string.len()+ 4 as usize, res_string.len());
    }
    #[test]
    fn test_correct_translation_6()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|5G>5GTEST|1936821C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|10V>10VECT|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(ref_string.len()+ 7 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_7()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|5G>5GTEST|1936821C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|10V>10VECT|1936821C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|14R>14RAPID|1936821C>T".to_string(),
        ];
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(ref_string.len()+ 11 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_8()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "frameshift|MAD1L1|ENST00000406869|protein_coding|-|10V>10VTESTFRAMESHIFT|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(24 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_9()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_deletion|MAD1L1|ENST00000406869|protein_coding|-|10VLSTLR>10V|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(33 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_10()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_deletion|MAD1L1|ENST00000406869|protein_coding|-|10VLSTLR>10R|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(33 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_11()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_deletion|MAD1L1|ENST00000406869|protein_coding|-|10VLSTLR>10R|1936821C>T".to_string(),
            "inframe_deletion|MAD1L1|ENST00000406869|protein_coding|-|28GSGLE>28E|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(29 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_12()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "stop_gained|MAD1L1|ENST00000406869|protein_coding|-|37G>37*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(36 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_13()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "stop_lost|MAD1L1|ENST00000406869|protein_coding|-|39*>39TEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(42 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_14()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "stop_lost|MAD1L1|ENST00000406869|protein_coding|-|38G*>39TEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(42 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_15()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*missense&inframe_altering|MAD1L1|ENST00000406869|protein_coding|-|34LERGG>34LTEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(38 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_16()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*frameshift&stop_retained|MAD1L1|ENST00000406869|protein_coding|-|20FISQRVEGGSGLEELERGG*>20LTEST*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(24 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_17()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*frameshift&stop_retained|MAD1L1|ENST00000406869|protein_coding|-|20FISQRVEGGSGLEELERGG*>20TEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(23 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_18()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "*stop_gained&inframe_altering|MAD1L1|ENST00000406869|protein_coding|-|20FISQRVEGGSGLEELERGG*>20|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(19 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_19()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "frameshift&stop_retained|MAD1L1|ENST00000406869|protein_coding|-|20FISQRVEGGSGLEELERGG*>20FLTESTTWO*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(28 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_20()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_deletion&stop_retained|MAD1L1|ENST00000406869|protein_coding|-|38*>38*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(38 as usize, res_string.len());
    }
    #[test]
    fn test_correct_translation_21()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "inframe_insertion&stop_retained|MAD1L1|ENST00000406869|protein_coding|-|38*>38*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(38 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_22()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "stop_gained&inframe_altering|MAD1L1|ENST00000406869|protein_coding|-|20FISQRVEGGSGLEELERGG*>20*|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(19 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_23()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "stop_lost&frameshift|MAD1L1|ENST00000406869|protein_coding|-|39*>39TEST|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(42 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_24()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "missense&inframe_altering|MAD1L1|ENST00000406869|protein_coding|-|34ERGG>34YEAP|1936821C>T".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLEELERGG".to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",res_string);
        assert_eq!(38 as usize, res_string.len());
    } 
    #[test]
    fn test_correct_translation_25()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "missense|MAD1L1|ENST00000265854|protein_coding|-|710E>710K|1816099C>T".to_string(),
            "missense|MAD1L1|ENST00000399654|protein_coding|-|706S>706L|1816110G>A".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|696R>696L|1816140C>A".to_string(),
            "missense|MAD1L1|ENST00000406869|protein_coding|-|666K>666N|1898200C>G".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|650R>650Q|1898249C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000455998|protein_coding|-|143L>143KQEL|2219360G>GCTCCTGCTT".to_string(),
            //"inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|190L>190KQEL|2219360G>GCTCCTGCT".to_string(),
            //"inframe_insertion|MAD1L1|ENST00000402746|protein_coding|-|98L>98KQEL|2219360G>GCTCCTGCTT".to_string(),
            //"inframe_deletion|MAD1L1|ENST00000399654|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string(),
            //"inframe_deletion|MAD1L1|ENST00000265854|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string()
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(721 as usize, res_string.len());
        assert_eq!(res_array[712],'K');
        assert_eq!(res_array[708],'L');
        assert_eq!(res_array[698],'L');
        assert_eq!(res_array[668],'N');
        assert_eq!(res_array[652],'Q');
    }  
    #[test]
    fn test_correct_translation_26()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "missense|MAD1L1|ENST00000265854|protein_coding|-|710E>710K|1816099C>T".to_string(),
            "missense|MAD1L1|ENST00000399654|protein_coding|-|706S>706L|1816110G>A".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|696R>696L|1816140C>A".to_string(),
            "missense|MAD1L1|ENST00000406869|protein_coding|-|666K>666N|1898200C>G".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|650R>650Q|1898249C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000455998|protein_coding|-|143L>143KQEL|2219360G>GCTCCTGCTT".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|190L>190KQEL|2219360G>GCTCCTGCT".to_string(),
            //"inframe_deletion|MAD1L1|ENST00000399654|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string(),
            //"inframe_deletion|MAD1L1|ENST00000265854|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string()
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(724 as usize, res_string.len());
        assert_eq!(res_array[715],'K');
        assert_eq!(res_array[711],'L');
        assert_eq!(res_array[701],'L');
        assert_eq!(res_array[671],'N');
        assert_eq!(res_array[655],'Q');
    }
    #[test]
    fn test_correct_translation_27()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "missense|MAD1L1|ENST00000265854|protein_coding|-|710E>710K|1816099C>T".to_string(),
            "missense|MAD1L1|ENST00000399654|protein_coding|-|706S>706L|1816110G>A".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|696R>696L|1816140C>A".to_string(),
            "missense|MAD1L1|ENST00000406869|protein_coding|-|666K>666N|1898200C>G".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|650R>650Q|1898249C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000455998|protein_coding|-|143L>143KQEL|2219360G>GCTCCTGCTT".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|190L>190KQEL|2219360G>GCTCCTGCT".to_string(),
            "inframe_deletion|MAD1L1|ENST00000399654|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string(),
            //"inframe_deletion|MAD1L1|ENST00000265854|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string()
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(722 as usize, res_string.len());
        assert_eq!(res_array[713],'K');
        assert_eq!(res_array[709],'L');
        assert_eq!(res_array[699],'L');
        assert_eq!(res_array[669],'N');
        assert_eq!(res_array[655],'Q');
    } 
    #[test]
    fn test_correct_translation_28()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
            "missense|MAD1L1|ENST00000265854|protein_coding|-|710E>710K|1816099C>T".to_string(),
            "missense|MAD1L1|ENST00000399654|protein_coding|-|706S>706L|1816110G>A".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|696R>696L|1816140C>A".to_string(),
            "missense|MAD1L1|ENST00000406869|protein_coding|-|666K>666N|1898200C>G".to_string(),
            "missense|MAD1L1|ENST00000265854|protein_coding|-|650R>650Q|1898249C>T".to_string(),
            "inframe_insertion|MAD1L1|ENST00000455998|protein_coding|-|143L>143KQEL|2219360G>GCTCCTGCTT".to_string(),
            "inframe_insertion|MAD1L1|ENST00000406869|protein_coding|-|190L>190KQEL|2219360G>GCTCCTGCT".to_string(),
            "inframe_deletion|MAD1L1|ENST00000437877|protein_coding|-|117DCL>117L|1898211GGCAGTC>G".to_string(),
            "inframe_deletion|MAD1L1|ENST00000399654|protein_coding|-|661DCL>661L|1898211GGCAGTC>G".to_string(),
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(720 as usize, res_string.len());
        assert_eq!(res_array[711],'K');
        assert_eq!(res_array[707],'L');
        assert_eq!(res_array[697],'L');
        assert_eq!(res_array[667],'N');
        assert_eq!(res_array[653],'Q');
    }       
    #[test]
    fn test_correct_translation_29()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
           "frameshift|MAD1L1|ENST00000406869|protein_coding|-|319RLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA*>319GETGPDHGPEHQDSRRPFQIRG*|1936821C>T+2213243T>TCTCC".to_string()
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(340 as usize, res_string.len());
    }       
    #[test]
    fn test_correct_translation_30()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec![
           "stop_gained|MAD1L1|ENST00000406869|protein_coding|-|82R>82*|2225457G>A".to_string()
        ]; 
        let alt_transcript= vcf_ds::AltTranscript::new(name, mutations);
        println!("{:#?}",alt_transcript); 
        let mut reference=HashMap::new(); 
        let ref_seq_array="MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA";
        reference.insert("ENST00000406869".to_string(),"MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQLEERAEQIRSKSHLIQVEREKMQMELSHKRARVELERAASTSARNYEREVDRNQELLTRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQASQEARADHEQQIKDLEQKLSLQEQDAAIVKNMKSELVRLPRLERELKQLREESAHLREMRETNGLLQEELEGLQRKLGRQEKMQETLVGLELENERLLAKLQSWERLDQTMGLSIRTPEDLSRFVVELQQRELALKDKNSAVTSSARGLEKARQQLQEELRQVSGQLLEERKKRETHEALARRLQKRVLLLTKERDGMRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAEMEAQLSQALEELGGQKQRADMLEMELKMLKSQSSSAEQSFLFSREEADTLRLKVEELEGERSRLEEEKRMLEAQLERRALQGDYDQSRTKVLHMSLNPTSVARQRLREDHSQLQAECERLRGLLRAMERGGTVPADLEAAAASLPSSKEVAELKKQVESAELKNQRLKEVFQTKIQEFRKACYTLTGYQIDITTENQYRLTSLYAEHPGDCLIFKATSPSGSKMQLLETEFSHTVGELIEVHLRRQDSIPAFLSSLTLELFSRQTVA".to_string());
        let res=TranscriptInstruction::from_alt_transcript(alt_transcript, &reference).unwrap(); 
        let test_gir=res.get_g_rep(&reference); 
        println!("{:#?}",test_gir); 
        let (res_array, _)=test_gir.execute(Engine::ST);
        let ref_string=ref_seq_array.to_string();
        println!("Input Sequence is:  ==>{:#?}",&ref_string);
        let res_string=res_array.iter().collect::<String>();
        println!("Result sequence is: ==>{:#?}",&res_string);
        assert_eq!(81 as usize, res_string.len());
    }       
}