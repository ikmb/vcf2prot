use std::collections::HashMap;
/// ## Definition 
/// The class act as a convient API for handling Fasta files, internally it utilizes a hashmap to associate 
/// every sequence id with it's sequence. 
/// ## Example 
///``` 
/// use ppgg_rust::data_structures::FastaFile; 
/// use ppgg_rust::readers::read_fasta_file; 
/// use std::path::Path; 
/// let path2file=Path::new("test_data/test_fasta_data1.fasta");
/// let fasta_file=read_fasta_file(path2file).unwrap(); 
/// assert_eq!(fasta_file.get_records().len(),3); 
/// assert!(fasta_file.is_in_records(&"seq1".to_string()));
///``` 
#[derive(Debug,Clone)]
pub struct FastaFile
{
    fastarecords:HashMap<String,String>
}
impl FastaFile
{
    /// ## Definition 
    /// Create a new instance from an input HashMap 
    /// ## Example 
    ///``` 
    /// use std::collections::HashMap; 
    /// use ppgg_rust::data_structures::FastaFile; 
    /// let mut test_map=HashMap::new(); 
    /// test_map.insert("seq1".to_string(),"test_seq".to_string()); 
    /// let fasta_file=FastaFile::FastaFile::new(test_map);
    /// assert!(fasta_file.is_in_records(&"seq1".to_string()));
    ///``` 
    pub fn new(fastarecords:HashMap<String,String>)->Self
    {
        FastaFile{fastarecords}
    }
    /// ## Definition 
    /// return a read only reference to the sequecne of the provided sequecne ID 
    /// ## Example 
    ///``` 
    /// use std::collections::HashMap; 
    /// use ppgg_rust::data_structures::FastaFile;
    /// let mut test_map=HashMap::new(); 
    /// test_map.insert("seq1".to_string(),"test_seq".to_string()); 
    /// let fasta_file=FastaFile::FastaFile::new(test_map);
    /// println!("The sequence of seq1 is  {}",fasta_file.get_record(&"seq1".to_string()).unwrap()); 
    ///``` 
    pub fn get_record(&self, seq_name:&String)->Result<&String,String>
    {
        if self.is_in_records(seq_name)
        {
            return Ok(&self.fastarecords.get(seq_name).unwrap()); 
        }
        Err(format!("The provided sequence name: {} not in not defined in the current fasta file",seq_name))
    }
    /// ## Definition 
    /// Return a read only reference to the file hashmap 
    /// ## Example 
    ///``` 
    /// use std::collections::HashMap; 
    /// use ppgg_rust::data_structures::FastaFile;
    /// let mut test_map=HashMap::new(); 
    /// test_map.insert("seq1".to_string(),"test_seq".to_string()); 
    /// let fasta_file=FastaFile::FastaFile::new(test_map);
    /// println!("The content of the Fasta file is {:#?}",fasta_file.get_records()); 
    ///``` 
    pub fn get_records(&self)->&HashMap<String,String>
    {
        &self.fastarecords
    }
    /// ## Definition 
    /// An indicator function, return True if the proivded seuqnece is defined inside the hash function and False otherwise.
    /// ## Example 
    ///``` 
    /// use std::collections::HashMap; 
    /// use ppgg_rust::data_structures::FastaFile;
    /// let mut test_map=HashMap::new(); 
    /// test_map.insert("seq1".to_string(),"test_seq".to_string()); 
    /// let fasta_file=FastaFile::FastaFile::new(test_map);
    /// assert!(fasta_file.is_in_records(&"seq1".to_string()));
    ///``` 
    pub fn is_in_records(&self, seq_name:&String)->bool
    {
        match self.fastarecords.get(seq_name)
        {
            Some(_)=>true,
            None=>false
        }
    }
    pub fn consume_and_get_hash_map(self)->HashMap<String,String>
    {
        self.fastarecords
    }
}