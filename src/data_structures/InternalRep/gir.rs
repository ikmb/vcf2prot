use std::collections::HashMap; 
use super::task::Task; 
use super::engines::Engine; 

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