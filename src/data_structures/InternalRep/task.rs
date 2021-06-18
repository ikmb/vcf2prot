// use a caret to load the data 

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct Task
{
    exe_code:u8,
    start_pos:usize,
    length:usize,
    start_pos_res:usize,
}
impl Task 
{
    /// Construct a new task 
    pub fn new(exe_code:u8, start_pos:usize,length:usize,
        start_pos_res:usize )->Self
    {
        Task{exe_code,start_pos,length,start_pos_res}
    }
    /// Execute the task of the two input streams ans the resulting vector of chars 
    /// ## Example  
    ///```     
    /// use ppgg_rust::data_structures::InternalRep::task::Task; 
    /// let test_stream_ref="ABCFEFGH"
    /// .chars()
    /// .collect::<Vec<char>>(); 
    /// let test_stream_alt=test_stream_ref.iter()
    ///    .rev()
    ///    .map(|c|c.clone())
    ///    .collect::<Vec<char>>(); 
    ///    let mut test_results=vec!['x';10];
    /// let mut expected_res=vec!['x';10];
    /// expected_res[8]='B'; 
    /// let task=Task::new(0,1,1,8); 
    /// task.execute(&mut test_results, &test_stream_ref, &test_stream_alt);
    /// assert_eq!(*test_results,*expected_res);
    ///``` 
    pub fn execute(&self, results_tape:&mut Vec<char>, ref_tape:&mut Vec<char>, alt_tape:&mut Vec<char>)
    {
        let end_bound_res=self.start_pos_res+self.length;
        let end_bound_stream=self.start_pos+self.length; 
        if self.exe_code==0
        {
            results_tape[self.start_pos_res..end_bound_res].clone_from_slice(&ref_tape[self.start_pos..end_bound_stream]); 
        }
        else
        {
            results_tape[self.start_pos_res..end_bound_res].clone_from_slice(&alt_tape[self.start_pos..end_bound_stream]); 
        }
    }
    pub fn get_mut_start_pos(&mut self)->&mut usize
    {
        &mut self.start_pos
    }
    pub fn get_mut_length(&mut self)->&mut usize
    {
        &mut self.length
    }
    pub fn get_length(&self)->usize
    {
        self.length
    }
    pub fn get_start_pos_res(&self)->usize
    {
        self.start_pos_res
    }
    pub fn get_mut_start_pos_res(&mut self)->&mut usize
    {
        &mut self.start_pos_res
    }
    pub fn get_execution_stream(&self)->&u8
    {
        &self.exe_code
    }
    pub fn shift_start_pos_stream(&mut self, num:&usize)
    {
        self.start_pos+=num; 
    }
    pub fn shift_start_pos_res(&mut self, num:&usize)
    {
        self.start_pos_res+=num; 
    }


}
#[cfg(test)]
pub mod test_task
{
    use super::*; 
    #[test]
    fn test_execute()
    {
        let mut test_stream_ref="ABCFEFGH"
                            .chars()
                            .collect::<Vec<char>>(); 
        let mut test_stream_alt=test_stream_ref.iter()
                                .rev()
                                .map(|c|c.clone())
                                .collect::<Vec<char>>(); 
        let mut test_results=vec!['x';10];
        // define the input streams 
        println!("{:#?}",test_stream_ref);
        println!("{:#?}",test_stream_alt);
        println!("{:#?}",test_results);
        let task=Task::new(0,1,1,8); 
        task.execute(&mut test_results, &mut test_stream_ref, &mut test_stream_alt);
        println!("{:#?}",test_results);
        let mut expected_res=vec!['x';10];
        expected_res[8]='B'; 
        assert_eq!(*test_results,*expected_res);
        let task2=Task::new(0,4,1,4); 
        task2.execute(&mut test_results, &mut test_stream_ref, &mut test_stream_alt);
        expected_res[4]='E'; 
        assert_eq!(*test_results,*expected_res);
        let task3=Task::new(0,6,2,6); 
        task3.execute(&mut test_results, &mut test_stream_ref, &mut test_stream_alt);
        expected_res[6]='G'; 
        expected_res[7]='H'; 
        assert_eq!(*test_results,*expected_res);
        println!("{:#?}",expected_res);
        println!("{:#?}",test_results);
    }
}