use std::collections::HashMap;
use std::fs; 
use std::io::Write;
use std::path::Path; 
use super::engines::Engine;
use super::proband_instructions::ProbandInstruction;
use super::sequence_tape::SequenceTape; 
/// an abstraction for a personalized proteome, it contains the proband_name and the sequence tap which contain the mutated_sequences
#[derive(Debug,Clone)]
pub struct PersonalizedGenome
{
    proband_name:String,
    seq_tape1:SequenceTape,
    seq_tape2:SequenceTape,
}
impl PersonalizedGenome
{
    /// Create a new instance from a sequence tape and a proband name 
    pub fn new(proband_name:String,seq_tape1:SequenceTape,seq_tape2:SequenceTape)->Self
    {
        PersonalizedGenome{proband_name,seq_tape1,seq_tape2}
    }
    /// write the personlized proteome to the results directory 
    /// ## Example 
    ///``` 
    /// use ppgg_rust::data_structures::InternalRep::{sequence_tape::SequenceTape,personalized_genome::PersonalizedGenome}; 
    /// use std::collections::HashMap; 
    /// let code_string1="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
    /// let code_string2="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
    /// let proband_name="code".to_string();
    /// let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
    /// res_map.insert("1".to_string(), (0,4)); 
    /// res_map.insert("2".to_string(), (5,9)); 
    /// res_map.insert("3".to_string(), (10,14)); 
    /// let seq_tape1=SequenceTape::new(code_string1, res_map.clone()).unwrap(); // this panic incase of length mismatch 
    /// let seq_tape2=SequenceTape::new(code_string2, res_map).unwrap(); 
    /// let personalized_proteome=PersonalizedGenome::new(proband_name, seq_tape1, seq_tape2); 
    /// personalized_proteome.write("test_data".to_string()).unwrap()
    ///```     
    pub fn write(&self, outdir:&String)->Result<(),String>
    {
        let res_string=format!("{}/{}.fasta",outdir,self.proband_name);
        let res_path=Path::new(&res_string); 
        let mut file_handle=match fs::File::create(res_path)
        {
            Ok(file)=>file,
            Err(err_msg)=>return Err(format!("Could not create {} because {}",res_path.display(),err_msg))
        }; 
        // write the content of the first sequence tape
        for (key,_) in self.seq_tape1.get_annotation().iter()
        {
            write!(&mut file_handle,">{}\n{}\n", key, self.seq_tape1.get_seq(key).unwrap()).unwrap();
        }
        // write the content of the first sequence tape
        for (key,_) in self.seq_tape2.get_annotation().iter()
        {
            write!(&mut file_handle,">{}\n{}\n", key, self.seq_tape2.get_seq(key).unwrap()).unwrap();
        }
        Ok(())
    }
    /// ## Summary
    /// create a new summary from a proband instruction, a reference proteome and an execution engine 
    pub fn from_proband_instruction(mut proband_instruction:ProbandInstruction, engine:Engine, ref_seq:&HashMap<String,String>)->Self
    {
        let proband_name=proband_instruction.proband_name; 
        let (res_1,annotations1)=proband_instruction.haplotype1_instruction.get_g_rep(ref_seq, engine.clone()).execute(engine.clone()); 
        let (res_2,annotations2)=proband_instruction.haplotype2_instruction.get_g_rep(ref_seq, engine.clone()).execute(engine.clone());
        let seq_tape1=SequenceTape::new(res_1.iter().collect::<String>(), annotations1).unwrap(); 
        let seq_tape2=SequenceTape::new(res_2.iter().collect::<String>(), annotations2).unwrap();
        PersonalizedGenome::new(proband_name, seq_tape1, seq_tape2) 
    }
}
#[cfg(test)]
mod test_personalized_proteome
{
    use super::*; 
    use std::collections::HashMap;
    #[test]
    pub fn test_personalized_proteome()->Result<(),String>
    {
        let code_string1="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
        let code_string2="SEQ1_SEQ2_SEQ3_SEQ4_SEQ5_SEQ6".to_string(); 
        let proband_name="code".to_string();
        let mut res_map:HashMap<String,(usize,usize)>=HashMap::new();
        res_map.insert("1".to_string(), (0,4)); 
        res_map.insert("2".to_string(), (5,9)); 
        res_map.insert("3".to_string(), (10,14)); 
        let seq_tape1=SequenceTape::new(code_string1, res_map.clone()).unwrap(); // this panic incase of length mismatch 
        let seq_tape2=SequenceTape::new(code_string2, res_map).unwrap(); 
        let personalized_proteome=PersonalizedGenome::new(proband_name, seq_tape1, seq_tape2); 
        personalized_proteome.write(&"test_data".to_string())
    }
}