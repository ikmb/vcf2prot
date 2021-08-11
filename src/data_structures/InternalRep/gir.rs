// load the modules and crate library 
use std::collections::HashMap; 
use super::task::Task; 
use super::engines::Engine; 
use crate::binders::binderCUDA; 


/// GIRL: Genomic intermediate representation language (GIRL) which us derived from sequence intermediate representation (SIR)
/// a generic representation for sequence editing tasks, it is composite of 
/// 1- a vector of tasks, see(Tasks for more details)
/// 2- an hashmap containing the bodundires of the resulting hashmap array 
/// 3- alt_stream: a vector of chars containing alterations, i.e. mutated amino acids 
/// 4- ref_stream: a vector of chars containing the reference stream
/// 5- res_array: a vector of chars containing the resulting arrays
/// the struct derives the Debug and the clone traits 
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
    /// ## Summary 
    /// Create a new instance using the main 5 units needed for representation 
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation:HashMap<String,(usize,usize)>=HashMap::new(); 
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// println!("The generated GIRL is: {:#?}",res);
    ///```
    pub fn new(g_rep:Vec<Task>, annotation:HashMap<String,(usize,usize)>, 
            alt_stream:Vec<char>, ref_stream:Vec<char>, res_array:Vec<char> )->Self
    {
        GIR{g_rep,annotation,alt_stream,ref_stream,res_array}
    }
    /// ## Summary
    /// Return a reference to the instance vector of tasks
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation:HashMap<String,(usize,usize)>=HashMap::new(); 
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // return a reference to the instance vector of tasks 
    /// let task = res.get_tasks();
    /// println!("The instance's tasks are: {:?}",task);
    ///´´´
    pub fn get_tasks(&self)->&Vec<Task>
    {
        &self.g_rep
    }
    /// ## Summary 
    /// Consume the instance and move its vector of tasks to the caller  
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation:HashMap<String,(usize,usize)>=HashMap::new(); 
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // return a reference to the instance vector of tasks 
    /// let task = res.consume_and_get_tasks();
    /// println!("The instance's tasks are: {:?}",task);
    /// //println!("This line will cause an error if printed as res has been consumed: {:#?}",res);
    ///´´´
    pub fn consume_and_get_tasks(self)->Vec<Task>
    {
        self.g_rep
    }
    /// ## Summary 
    /// Consume the instance and return its component as a tuple of elements  
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation:HashMap<String,(usize,usize)>=HashMap::new(); 
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // return a reference to the instance vector of tasks 
    /// let task = res.consume_and_get_tasks();
    /// println!("The instance's tasks are: {:?}",task);
    /// //println!("This line will cause an error if printed as res has been consumed: {:#?}",res);
    pub fn consumer_and_get_resources(self)->(Vec<Task>,HashMap<String,(usize,usize)>,Vec<char>,Vec<char>,Vec<char>)
    {
        (self.g_rep,self.annotation,self.alt_stream,self.ref_stream,self.res_array)
    }
    /// ## Summary 
    /// Return a hash map containing the annotations associated in the reference array 
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation:HashMap<String,(usize,usize)>=HashMap::new(); 
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // return a reference to the instance vector of tasks 
    /// let annotations = res.get_annotation();
    /// println!("The instance's annotation is: {:?}",annotations);
    ///```
    pub fn get_annotation(&self)->&HashMap<String,(usize,usize)>
    {
        &self.annotation
    }
    /// ## Summary 
    /// Return the max index in the tasks vector
    /// ## Examples 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::task; 
    /// use std::collections::HashMap; 
    /// // let's define some place holders for a task 
    /// let g_rep:Vec<Task>=Vec::new(); 
    /// let annotation=HashMap::new(); 
    /// // let's add some annotations to the code 
    /// annotations.insert("Sequence_1".to_string(), (0 as usize,25 as usize));
    /// annotations.insert("Sequence_2".to_string(), (25 as usize, 50 as usize))
    /// let alt_stream:Vec<char>=Vec::new(); 
    /// let ref_stream:Vec<char>=Vec::new(); 
    /// let res_array:Vec<char> =Vec::new(); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // return a reference to the instance vector of tasks 
    /// assert_eq!(res.get_results_max(),50 as usize);
    ///``` 
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
    /// ## Summary 
    /// execute and consume the representation to return a vector of chars containing the edited sequences 
    /// along with a hashmap containing index of features in the results vector 
    /// ## Example 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::{task,engines}; 
    /// use std::collections::HashMap;
    /// // let's define some dummy example data
    /// let g_rep=Vec::new(); 
    /// let g_rep.push(task::Task::new(0/* execution_code is zero, i.e. reference stream */,
    ///                 0 /* start copying from position 0 */,
    ///                 4 /* copies 4 amino acids*/, 
    ///                 0 /* insert at position 0 in the reference array */)); 
    /// let g_rep.push(task::Task::new(1 /* execution_code is 1, i.e. alternative stream */,
    ///                 0 /* start copying from position 0 */,
    ///                 1 /* copies 1 amino acids*/, 
    ///                 4 /* insert at position 0 in the reference array */)); 
    /// let annotation=HashMap::new(); 
    /// annotation.insert("Seq_1".to_string(),(0 as usize,5 as usize)); 
    /// let alt_stream=vec!['G']; 
    /// let ref_stream=vec!['T','E','S','T']; 
    /// let res_array:Vec<char> =Vec::with_capacity(5); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // execute the GIR with a single threaded engine 
    /// let (res_char_array, res_hashmap)=res.execute(engines::Engine::from_str("st")); 
    /// println!("Results array: {:#?}",res_char_array); 
    /// println!("Result hashmap is: {:#?}", res_hashmap);
    ///``` 
    pub fn execute(self, engine:Engine)->(Vec<char>,HashMap<String,(usize,usize)>)
    {        
        match engine 
        {
            Engine::ST | Engine::MT =>
            {
                match std::env::var("DEBUG_CPU_EXEC")
                {
                    Ok(_)=>
                    {
                        println!("Validating the execution tasks on the CPU engine ....");
                        for idx in 1..self.g_rep.len()
                        {
                            if self.g_rep[idx].get_start_pos_res()!=self.g_rep[idx-1].get_start_pos_res() + self.g_rep[idx-1].get_length()
                            {
                                println!("************ CPU Execution Table *********");
                                println!("index\tstream\tstart_position\tlength\tposition_results\t");
                                for idx in 0..self.g_rep.len()
                                {
                                    println!(
                                        "{}\t{}\t{}\t{}\t{}\t",idx,self.g_rep[idx].get_execution_stream(),
                                        self.g_rep[idx].get_start_pos(),
                                        self.g_rep[idx].get_length(),
                                        self.g_rep[idx].get_start_pos_res()
                                        );
                                }
                                panic!("Critical failure in the calculations was encountered: position: {} the sum {} does not equal previous inputs: {} and {} \n",
                                idx,self.g_rep[idx].get_start_pos_res(),self.g_rep[idx-1].get_start_pos_res(),self.g_rep[idx-1].get_length());
                            }
                        }
                    },
                    Err(_)=>()
                }
                let mut res_array=self.res_array; 
                let mut ref_stream=self.ref_stream;
                let mut alt_stream=self.alt_stream;
                self.g_rep.iter().for_each(|task| task.execute(&mut res_array, &mut ref_stream, &mut alt_stream));
                (res_array,self.annotation)
            },
            Engine::GPU => 
            {
                let (exec_code,start_pos,length,start_pos_res,res_array,
                    ref_array,alt_array,  annotation )= self.consume_and_produce_produce_content(); 
                // cast as u8; define the results array 
                let mut res_array=res_array.into_iter().map(|val|val as u8).collect::<Vec<_>>(); 
                let ref_array=ref_array.into_iter().map(|val|val as u8).collect::<Vec<_>>(); 
                let alt_array=alt_array.into_iter().map(|val|val as u8).collect::<Vec<_>>(); 
                let err_code; 
                // validate the execution before  
                match std::env::var("DEBUG_GPU")
                {
                    Ok(_)=>
                    {
                        println!("Validating the execution tasks ....");
                        for idx in 1..exec_code.len()
                        {
                            if start_pos_res[idx]!=start_pos_res[idx-1]+length[idx-1]
                            {
                                println!("************ GPU Execution Table *********");
                                println!("index\tstream\tstart_position\tlength\tposition_results\t");
                                for idx in 0..exec_code.len()
                                {
                                    println!("{}\t{}\t{}\t{}\t{}\t",idx,exec_code[idx],start_pos[idx],length[idx],start_pos_res[idx]);
                                }
                                panic!("Critical failure in the calculations was encountered: position: {} the sum {} does not equal previous inputs: {} and {} \n",
                                idx,start_pos_res[idx],start_pos_res[idx-1],length[idx-1]);
                            }
                        }
                    },
                    Err(_)=>()
                }
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
                let res_array=res_array.into_iter().map(|elem| elem as char).collect::<Vec<_>>(); 
                (res_array, annotation)
            }
        }
    }   
    /// ## Summary  ,ref_array,alt_array,annotation)
    /// Consume the instance and return the following arrays:
    /// 1. A vector of usize containing the execution code 
    /// 2. A vector of usize containing the start positions in reference
    /// 3. A vector of usize containing the length of each instruction  
    /// 4. A char array containing the results array
    /// 5. A char array containing the reference array 
    /// 6. A char array containing the alternative array 
    /// 7. A hashmap array containing a map between a string and a tuple of two usize integers 
    /// ## Example 
    /// ```rust
    /// // let's load the need modules first 
    /// use ppgg_rust::data_structures::internal_representation::{task,engines}; 
    /// use std::collections::HashMap;
    /// // let's define some dummy example data
    /// let g_rep=Vec::new(); 
    /// let g_rep.push(task::Task::new(0/* execution_code is zero, i.e. reference stream */,
    ///                 0 /* start copying from position 0 */,
    ///                 4 /* copies 4 amino acids*/, 
    ///                 0 /* insert at position 0 in the reference array */)); 
    /// let g_rep.push(task::Task::new(1 /* execution_code is 1, i.e. alternative stream */,
    ///                 0 /* start copying from position 0 */,
    ///                 1 /* copies 1 amino acids*/, 
    ///                 4 /* insert at position 0 in the reference array */)); 
    /// let annotation=HashMap::new(); 
    /// annotation.insert("Seq_1".to_string(),(0 as usize,5 as usize)); 
    /// let alt_stream=vec!['G']; 
    /// let ref_stream=vec!['T','E','S','T']; 
    /// let res_array:Vec<char> =Vec::with_capacity(5); 
    /// let res=GIR::new(g_rep, annotation, alt_stream, ref_stream, res_array); 
    /// // execute the GIR with a single threaded engine 
    /// let (exec_code,start_pos,length,start_pos_res,res_array,ref_array,alt_array,annotation)=res.consume_and_produce_produce_content(); 
    /// println!("The execution code is {:#?}",exec_code); 
    /// println!("The execution code is {:#?}",start_pos);
    /// println!("The length is {:#?}",length); 
    /// println!("The start position in results is {:#?}",start_pos_res);  
    /// println!("The results array is {:#?}",res_array);  
    /// println!("The reference array is {:#?}",ref_array);  
    /// println!("The alternative array is {:#?}",alt_array);  
    /// println!("The annotation map is {:#?}",annotation);  
    ///``` 
    fn consume_and_produce_produce_content(self)->(Vec<usize>,Vec<usize>,Vec<usize>, 
        Vec<usize>, Vec<char>, Vec<char>, Vec<char>, HashMap<String,(usize,usize)>)
    {
        let mut exec_code=Vec::with_capacity(self.g_rep.len()); 
        let mut start_pos=Vec::with_capacity(self.g_rep.len()); 
        let mut length=Vec::with_capacity(self.g_rep.len()); 
        let mut start_pos_res=Vec::with_capacity(self.g_rep.len()); 
        for task in self.g_rep.iter()
        {
            exec_code.push(task.get_stream() as usize); 
            start_pos.push(task.get_start_pos()); 
            length.push(task.get_length()); 
            start_pos_res.push(task.get_start_pos_res())
        }
        let ( res_array,  ref_array,  alt_array, annotation)=(self.res_array, self.ref_stream, self.alt_stream, self.annotation); 
        (exec_code,start_pos,length,start_pos_res,res_array,ref_array,alt_array,annotation)
    }
}