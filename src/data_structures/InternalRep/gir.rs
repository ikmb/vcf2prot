use std::collections::HashMap; 
use super::task::Task; 
use super::engines::Engine; 
use rayon::prelude::*; 

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
    pub fn new(g_rep:Vec<Task>, annotation:HashMap<String,(usize,usize)>, 
            alt_stream:Vec<char>, ref_stream:Vec<char>, res_array:Vec<char> )->Self
    {
        GIR{g_rep,annotation,alt_stream,ref_stream,res_array}
    }
    pub fn add(&mut self, taks:Task, seq_name:String, seq_pos:(usize,usize))
    {
        self.g_rep.push(taks);
        self.annotation.insert(seq_name, seq_pos); 
    }
    pub fn get_taks(&self)->&Vec<Task>
    {
        &self.g_rep
    }
    pub fn consume_and_get_tasks(mut self)->Vec<Task>
    {
        self.g_rep
    }
    pub fn consumer_and_get_resources(mut self)->(Vec<Task>,HashMap<String,(usize,usize)>,Vec<char>,Vec<char>,Vec<char>)
    {
        (self.g_rep,self.annotation,self.alt_stream,self.ref_stream,self.res_array)
    }
    pub fn get_annotation(&self)->&HashMap<String,(usize,usize)>
    {
        &self.annotation
    }
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
    pub fn execute(self, engine:Engine)->(Vec<char>,HashMap<String,(usize,usize)>)
    {
        match engine 
        {
            Engine::ST=>
            {
                let mut res_array=self.res_array; 
                let mut ref_stream=self.ref_stream;
                let mut alt_stream=self.alt_stream;
                self.g_rep.iter().for_each(|task| task.execute(&mut res_array, &mut ref_stream, &mut alt_stream));
                (res_array,self.annotation)
                
            },
            _=>panic!("Not yet, not yet for parallel code")

        }
    }

}