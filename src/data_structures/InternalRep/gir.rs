use std::cmp;
use std::collections::HashMap; 
use std::mem::swap;
use std::sync::{Arc, Mutex};
use super::task::Task; 
use super::engines::Engine; 
use rayon::prelude::*; 
use scoped_threadpool;
use crossbeam; 
use num_cpus; 
use crate::binders::binderCUDA; 

/// GIRL: Genomic intermediate representation language (GIRL)
/// a generic representation for sequence editing tasks, it is composite of 
/// 1- a vector of tasks, see(Tasks for more details)
/// 2- an hashmap containing the bodundires of the resulting hashmap array 
/// 3- alt_stream: a vector of chars containing alterations, i.e. mutated amino acids 
/// 4- ref_stream: a vector of chars containing the reference stream
/// 5- res_array: a vector of chars containing the resulting arrays
#[derive(Debug,Clone)]
pub struct GIR
{
    g_rep:Vec<Task>,
    annotation:HashMap<String,(usize,usize)>, 
    alt_stream:Vec<char>,
    ref_stream:Vec<char>,
    res_array:Vec<char>
}
impl GIR
{
    /// create a new instance using the main 5 units needed for representation 
    pub fn new(g_rep:Vec<Task>, annotation:HashMap<String,(usize,usize)>, 
            alt_stream:Vec<char>, ref_stream:Vec<char>, res_array:Vec<char> )->Self
    {
        GIR{g_rep,annotation,alt_stream,ref_stream,res_array}
    }
    /// return a reference to the instance vector of tasks 
    pub fn get_taks(&self)->&Vec<Task>
    {
        &self.g_rep
    }
    /// consume the instance and return move its vector of tasks to the caller  
    pub fn consume_and_get_tasks(self)->Vec<Task>
    {
        self.g_rep
    }
    /// consume the instance and return its compoanant as a tuple of elements  
    pub fn consumer_and_get_resources(self)->(Vec<Task>,HashMap<String,(usize,usize)>,Vec<char>,Vec<char>,Vec<char>)
    {
        (self.g_rep,self.annotation,self.alt_stream,self.ref_stream,self.res_array)
    }
    /// return a hash map containing the annotations associated in the reference array  
    pub fn get_annotation(&self)->&HashMap<String,(usize,usize)>
    {
        &self.annotation
    }
    /// return the max index in the tasks vector 
    pub fn get_results_max(&self)->usize
    {
        let mut max_length=0; 
        for (_,value) in self.annotation.iter()
        {
            if value.1>max_length
            {
                max_length=value.1
            }
        }
        max_length
    }
    /// execute and consume the representation to return a vector of chars containing the edited sequences 
    /// along with a hashmap containing index of features in the results vector 
    pub fn execute(mut self, engine:Engine)->(Vec<char>,HashMap<String,(usize,usize)>)
    {
        match engine 
        {
            Engine::ST | Engine::MT =>
            {
                let mut res_array=self.res_array; 
                let mut ref_stream=self.ref_stream;
                let mut alt_stream=self.alt_stream;
                self.g_rep.iter().for_each(|task| task.execute(&mut res_array, &mut ref_stream, &mut alt_stream));
                (res_array,self.annotation)
            },
            Engine::GPU => 
            {
                let (exec_code,start_pos,length,start_pos_res,mut res_array,
                    ref_array,alt_array,  annotation )= self.consume_and_produce_produce_content(); 
                
                let exec_code =exec_code.into_iter().map(|elem| elem as usize).collect::<Vec<_>>(); 
                let err_code; 
                unsafe
                {
                    err_code=binderCUDA::kernel_wrapper(res_array.as_mut_ptr(),
                    ref_array.as_ptr(),alt_array.as_ptr(),
                    exec_code.as_ptr(), start_pos.as_ptr(), length.as_ptr(), 
                    start_pos_res.as_ptr(), exec_code.len(), res_array.len(), 
                    ref_array.len(), alt_array.len()); 
                }
                match err_code
                {
                    0=>(), 
                    1=>panic!("Allocating arrays on the GPU failed"),
                    2=>panic!("Failure with copying the data to the GPU"), 
                    3=>panic!("Launching the kernel failed"), 
                    4=>panic!("Kernel execution failed"),
                    5=>panic!("Copying the results array to the host failed"),
                    _=>panic!("Unknown error was encountered")
                }
                (res_array, annotation)
            }
        }
    }

    fn consume_and_produce_produce_content(mut self)->(Vec<u8>,Vec<usize>,Vec<usize>, 
        Vec<usize>, Vec<char>, Vec<char>, Vec<char>, HashMap<String,(usize,usize)>)
    {
        let mut exec_code=Vec::with_capacity(self.g_rep.len()); 
        let mut start_pos=Vec::with_capacity(self.g_rep.len()); 
        let mut length=Vec::with_capacity(self.g_rep.len()); 
        let mut start_pos_res=Vec::with_capacity(self.g_rep.len()); 
        for task in self.g_rep.iter()
        {
            exec_code.push(task.get_stream()); 
            start_pos.push(task.get_start_pos()); 
            length.push(task.get_length()); 
            start_pos_res.push(task.get_start_pos_res())
        }
        let ( res_array,  ref_array,  alt_array, annotation)=(self.res_array, self.res_array, self.alt_stream, self.annotation); 
        (exec_code,start_pos,length,start_pos_res,res_array,ref_array,alt_array,annotation)
    }






















}