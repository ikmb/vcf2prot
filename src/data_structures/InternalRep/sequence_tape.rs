use std::collections::HashMap; 
use std::path::Path; 
use std::fs; 
use std::io::Write;
/// An abstraction for a sequence tape, where more than one sequence are annotated in an head to tail fashion 
/// and a has map that stores the sequence name and the boundries, i.e. the start and the end point in the sequence
/// are stored. 
#[derive(Debug,Clone)]
pub struct SequenceTape
{
    seq_str:String,
    annotations:HashMap<String,(usize,usize)>
}

impl SequenceTape
{
    /// create a new sequence map from a seuqnece tape and an annotation hash map 
    /// ## Example 
    ///``` 
    /// use ppgg_rust::data_structures::InternalRep::sequence_tape::SequenceTape; 
    /// use std::collections::HashMap; 
    /// let code_string="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
    /// let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
    /// res_map.insert("1".to_string(), (0,4)); 
    /// res_map.insert("2".to_string(), (5,9)); 
    /// res_map.insert("3".to_string(), (10,14)); 
    /// let seq_tape=SequenceTape::new(code_string, res_map).unwrap(); // this panic incase of length mismatch 
    /// //check the correct mapping between the annotations
    /// assert_eq!("SEQ1",seq_tape.get_seq(&"1".to_string()).unwrap()); 
    /// assert_eq!("SEQ2",seq_tape.get_seq(&"2".to_string()).unwrap());
    /// assert_eq!("SEQ3",seq_tape.get_seq(&"3".to_string()).unwrap());
    ///``` 
    pub fn new(seq_str:String,annotations:HashMap<String,(usize,usize)>)->Result<Self,String>
    {
        let max_index=SequenceTape::get_max_index(&annotations); 
        if max_index>seq_str.len()
        {
            return Err(format!("Bad Tape Encountered, the provided maximum index is {} while tape length is {} ",max_index,seq_str.len())); 
        }
        Ok(SequenceTape{seq_str,annotations})
    }
    /// Write the sequence tap to a fasta file on disk 
    /// ## Example 
    ///``` 
    /// use std::path::Path;
    /// use ppgg_rust::data_structures::InternalRep::sequence_tape::SequenceTape; 
    /// use std::collections::HashMap; 
    /// let code_string="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
    /// let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
    /// res_map.insert("1".to_string(), (0,4)); 
    /// res_map.insert("2".to_string(), (5,9)); 
    /// res_map.insert("3".to_string(), (10,14)); 
    /// let seq_tape=SequenceTape::new(code_string, res_map).unwrap(); // this panic incase of length mismatch 
    /// seq_tape.write_to_fasta(Path::new("test_data/test_file.fasta")).unwrap();
    ///``` 
    pub fn write_to_fasta(&self,output_file_name:&Path)->Result<(),String>
    {
        let mut file_handle=match fs::File::create(output_file_name)
        {
            Ok(file)=>file,
            Err(err_msg)=>return Err(format!("Could not create {} because {}",output_file_name.display(),err_msg))
        }; 
        for (key,_) in self.annotations.iter()
        {
            write!(&mut file_handle,">{}\n{}\n", key, self.get_seq(key).unwrap()).unwrap();
        }
        Ok(())
    }
    /// return the hash map containing the annotation hash map 
    pub fn get_annotation(&self)->&HashMap<String,(usize,usize)>
    {
        &self.annotations
    }
    /// return the sequence corresponding to the provided sequence name 
    pub fn get_seq(&self,seq_name:&String)->Result<&str,String>
    {
        let res=match self.annotations.get(seq_name) {
            Some(res)=>res,
            None=>return Err(format!("The provided sequence name: {}, is not defined in the current tabe",seq_name))
        };
        let res_string=&self.seq_str[res.0..res.1];
        Ok(res_string)
    }
    /// return the sequence corresponding to the maximum index of the tape
    /// ## Example 
    ///``` 
    /// use std::collections::HashMap;  
    /// use ppgg_rust::data_structures::InternalRep::sequence_tape::SequenceTape; 
    /// let code_string="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
    /// let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
    /// res_map.insert("1".to_string(), (0,4)); 
    /// res_map.insert("2".to_string(), (5,9)); 
    /// res_map.insert("3".to_string(), (10,14)); 
    /// res_map.insert("4".to_string(), (15,19)); 
    /// res_map.insert("5".to_string(), (20,24)); 
    /// res_map.insert("6".to_string(), (25,29)); 
    /// let seq_tape=SequenceTape::new(code_string, res_map).unwrap(); 
    /// assert_eq!(SequenceTape::get_max_index(seq_tape.get_annotation()),29);
    ///``` 
    pub fn get_max_index(annotation:&HashMap<String,(usize,usize)>)->usize
    {
        let mut max=0; 
        for (_,value) in annotation.iter()
        {
            if value.1 >max
            {
                max=value.1
            }
        }
        max
    }
}
#[cfg(test)]
pub mod test_sequence_tape_module
{
    use std::usize;
    use super::*;
    #[test]
    pub fn test_sequence_tape()
    {
        let code_string="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
        let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
        res_map.insert("1".to_string(), (0,4)); 
        res_map.insert("2".to_string(), (5,9)); 
        res_map.insert("3".to_string(), (10,14)); 
        res_map.insert("4".to_string(), (15,19)); 
        res_map.insert("5".to_string(), (20,24)); 
        res_map.insert("6".to_string(), (25,29)); 
        let seq_tape=SequenceTape::new(code_string, res_map).unwrap(); 
        // check the correct mapping between the annotations
        assert_eq!("SEQ1",seq_tape.get_seq(&"1".to_string()).unwrap()); 
        assert_eq!("SEQ2",seq_tape.get_seq(&"2".to_string()).unwrap());
        assert_eq!("SEQ3",seq_tape.get_seq(&"3".to_string()).unwrap());
        assert_eq!("SEQ4",seq_tape.get_seq(&"4".to_string()).unwrap());
        assert_eq!("SEQ5",seq_tape.get_seq(&"5".to_string()).unwrap());
        assert_eq!("SEQ6",seq_tape.get_seq(&"6".to_string()).unwrap());
    }
    #[test]
    #[should_panic]
    pub fn test_sequence_tape2()
    {
        let code_string="SEQ1_SEQ2_SEQ3".to_string(); 
        let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
        res_map.insert("1".to_string(), (0,4)); 
        res_map.insert("2".to_string(), (5,9)); 
        res_map.insert("3".to_string(), (10,code_string.len()+1)); 
       
        SequenceTape::new(code_string, res_map).unwrap(); 
    }
    #[test]
    pub fn test_fasta_write()
    {
        let code_string="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
        let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
        res_map.insert("1".to_string(), (0,4)); 
        res_map.insert("2".to_string(), (5,9)); 
        res_map.insert("3".to_string(), (10,14)); 
        res_map.insert("4".to_string(), (15,19)); 
        res_map.insert("5".to_string(), (20,24)); 
        res_map.insert("6".to_string(), (25,29)); 
        let seq_tape=SequenceTape::new(code_string, res_map).unwrap(); 
        seq_tape.write_to_fasta(Path::new("test_data/test_file.fasta")).unwrap();
    }
    
}
