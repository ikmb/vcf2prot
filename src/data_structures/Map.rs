use std::collections::HashMap; 
use super::vcf_ds::AltTranscript; 
use super::InternalRep::proband_instructions::ProbandInstruction;
use serde::{Deserialize, Serialize};

#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct IntMap
{
    pub proband_name:String,
    mutations1:Vec<AltTranscript>,
    mutations2:Vec<AltTranscript>,
}
impl IntMap
{
    pub fn new(proband_name:String,mutations1:Vec<AltTranscript>,mutations2:Vec<AltTranscript>)->Self
    {
        IntMap{proband_name,mutations1,mutations2}
    }
    pub fn get_mutations_ref(&self)->(&Vec<AltTranscript>,&Vec<AltTranscript>)
    {
        (&self.mutations1,&self.mutations2)
    }
    pub fn get_name(&self)->&String
    {
        &self.proband_name
    } 
    pub fn consume_and_get_vecs(mut self)->(Vec<AltTranscript>,Vec<AltTranscript>)
    {
        (self.mutations1,self.mutations2)
    }
}

/// A data structure used to represent the early links between a map its mutations.
/// the struct owns three data strucutres: a proband_name which hold the name of the proband, 
/// mutations1 which holds all mutations in the first haplotype 
/// mutations2 which holds all mutations in the second haplotype 
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct EarlyMap
{
    proband_name:String,
    mutations1:Vec<String>,
    mutations2:Vec<String>,
}
impl EarlyMap
{
    /// create a new instance using a proband-name, a vector of string represent the consequences of the first haplotype 
    /// and a second vector represent the consequences of the second haplotypes 
    /// ## Examples 
    ///``` 
    /// use ppgg_rust::data_structures::Map::EarlyMap;
    /// let proband_name="Test_name".to_string(); 
    /// let mutations1="mutation1_1,mutation1_3".split(",").map(|elem| elem.to_string()).collect::<Vec<String>>(); 
    /// let mutations2="mutation1_2,mutation1_4".split(",").map(|elem| elem.to_string()).collect::<Vec<String>>(); 
    /// let results=EarlyMap::new(proband_name.clone(),mutations1,mutations2); 
    /// assert_eq!(*results.get_proband_name(),proband_name); 
    /// assert_eq!(results.get_mutations_ref(),(&"mutation1_1,mutation1_3".split(",").map(|elem| elem.to_string()).collect::<Vec<String>>(),
    /// &"mutation1_2,mutation1_4".split(",").map(|elem| elem.to_string()).collect::<Vec<String>>()))
    ///```
    pub fn new(proband_name:String, mutations1:Vec<String>, mutations2:Vec<String>)->Self
    {
        EarlyMap{proband_name,mutations1,mutations2}
    }
    /// create a new instance using a proband-name, and allocate two vectors to hold the generated strings, the expected number of mutations is 
    /// determined by the parameter, expected_number 
    /// ## Examples 
    ///``` 
    /// use ppgg_rust::data_structures::Map::EarlyMap;
    /// let proband_name="Test_name".to_string(); 
    /// let num_mutation=10; 
    /// let results=EarlyMap::with_capacity(proband_name.clone(),10); 
    /// assert_eq!(*results.get_proband_name(),proband_name); 
    /// assert_eq!(results.get_mutations_ref().0.capacity(),num_mutation);
    /// assert_eq!(results.get_mutations_ref().0.len(),0);
    /// assert_eq!(results.get_mutations_ref().1.capacity(),num_mutation);
    /// assert_eq!(results.get_mutations_ref().0.len(),0);
    ///```
    pub fn with_capacity(proband_name:String, expected_number:usize)->Self
    {
        let mut mutations1:Vec<String>= Vec::with_capacity(expected_number);
        let mut mutations2:Vec<String>= Vec::with_capacity(expected_number);
        EarlyMap{proband_name,mutations1,mutations2}
    }
    /// add a new mutation to the current instance
    /// ## Examples 
    ///``` 
    /// use ppgg_rust::data_structures::Map::EarlyMap;
    /// let proband_name="Test_name".to_string(); 
    /// let num_mutation=10; 
    /// let mut results=EarlyMap::with_capacity(proband_name.clone(),10); 
    /// results.add_mutation("mutation1".to_string(),1); 
    /// results.add_mutation("mutation2".to_string(),2);
    /// assert_eq!(*results.get_proband_name(),proband_name); 
    /// assert_eq!(results.get_mutations_ref(),(&vec!["mutation1".to_string()], &vec!["mutation2".to_string()]));
    ///```
    pub fn add_mutation(&mut self, mutation:String, haplotype:u8)
    {
        if haplotype==1
        {
            self.mutations1.push(mutation);
        }
        else 
        {
            self.mutations2.push(mutation);
        }
        
    }
    /// return a tuple that contain reference to the instances vector of string 
    /// ## Examples 
    ///``` 
    /// use ppgg_rust::data_structures::Map::EarlyMap;
    /// let proband_name="Test_name".to_string(); 
    /// let num_mutation=10; 
    /// let mut results=EarlyMap::with_capacity(proband_name.clone(),10); 
    /// results.add_mutation("mutation1".to_string(),1); 
    /// results.add_mutation("mutation2".to_string(),2);
    /// assert_eq!(*results.get_proband_name(),proband_name); 
    /// assert_eq!(results.get_mutations_ref(),(&vec!["mutation1".to_string()], &vec!["mutation2".to_string()]));
    ///```
    pub fn get_mutations_ref(&self)->(&Vec<String>,&Vec<String>)
    {
        (&self.mutations1,&self.mutations2)
    }
    /// return a reference to the proband names 
    /// ## Examples 
    ///``` 
    /// use ppgg_rust::data_structures::Map::EarlyMap;
    /// let proband_name="Test_name".to_string(); 
    /// let num_mutation=10; 
    /// let mut results=EarlyMap::with_capacity(proband_name.clone(),10); 
    /// assert_eq!(*results.get_proband_name(),proband_name); 
    ///```
    pub fn get_proband_name(&self)->&String
    {
        &self.proband_name
    }
}