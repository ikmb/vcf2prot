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
    /// ## Summary 
    /// Construct a new task 
    pub fn new(exe_code:u8, start_pos:usize,length:usize,
        start_pos_res:usize )->Self
    {
        Task{exe_code,start_pos,length,start_pos_res}
    }
    /// ## Summary
    /// Execute the task of the two input streams ans the resulting vector of chars 
    /// ## Example  
    ///```     
    /// use ppg_rust::data_structures::InternalRep::task::Task; 
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
    /// ## Summary
    /// get a mutable reference to the start position 
    pub fn get_mut_start_pos(&mut self)->&mut usize
    {
        &mut self.start_pos
    }
    /// ## Summary
    /// get a mutable reference to the instance's length 
    pub fn get_mut_length(&mut self)->&mut usize
    {
        &mut self.length
    }
    /// ## Summary
    ///  return the instance's length 
    pub fn get_length(&self)->usize
    {
        self.length
    }
    /// ## Summary
    ///  return the instance's start pos in the results array 
    pub fn get_start_pos_res(&self)->usize
    {
        self.start_pos_res
    }
    /// ## Summary
    ///  return a mutable reference to the start position in the results array 
    pub fn get_mut_start_pos_res(&mut self)->&mut usize
    {
        &mut self.start_pos_res
    }
    /// ## Summary
    ///  return the execution stream 
    pub fn get_execution_stream(&self)->&u8
    {
        &self.exe_code
    }
    /// ## Summary
    ///  return the execution stream 
    #[inline]
    pub fn get_stream(&self)->u8
    {
        self.exe_code
    }
    /// ## Summary
    ///  return the start position in the input stream 
    #[inline]
    pub fn get_start_pos(self)->usize
    {
        self.start_pos
    }
    /// ## Summary
    ///  shirt, i.e. change the start position in the stream 
    pub fn shift_start_pos_stream(&mut self, num:&usize)
    {
        self.start_pos+=*num; 
    }
    /// ## Summary
    ///  shirt, i.e. change the start position in the result array  
    pub fn shift_start_pos_res(&mut self, num:&usize)
    {
        self.start_pos_res+=*num; 
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
        let task=Task::new(0,1,1,8); 
        task.execute(&mut test_results, &mut test_stream_ref, &mut test_stream_alt);
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
    }
}