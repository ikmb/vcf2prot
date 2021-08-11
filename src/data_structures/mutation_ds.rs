use std::{cmp::Ordering, str::FromStr}; 
use crate::functions::text_parser; 
/// an enumarator that contain the supported mutation, namely, MisSense for missense mutations, 
/// InframeInsertion, i.e. inserions,  InframeDeletion, i.e deletion, FrameShift for frameshifts,
/// StopGain, i.e. stop gain and StopLost, i.e. stop lost.
/// it derives the Debug, Clone, PartialEq trait and impelement FromStr trait.
/// ## Example  
///``` 
/// let test_case="missense"; 
/// use std::str::FromStr;
/// use ppg_rust::data_structures::mutation_ds::MutationType;
/// assert_eq!(MutationType::MisSense, MutationType::from_str(test_case).unwrap());
///``` 
#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]
pub enum MutationType{
                    MisSense,InframeInsertion,InframeDeletion,FrameShift,StopGained,StopLost,
                    SMisSense,SInframeInsertion,SInframeDeletion,SFrameShift,SMisSenseAndInframeAltering,
                    SFrameShiftAndStopRetained,SStopGainedAndInframeAltering,FrameShiftAndStopRetained,
                    InframeDeletionAndStopRetained,InframeInsertionAndStopRetained,StopGainedAndInframeAltering,
                    StartLost,SStopGained,StopLostAndFrameShift, MissenseAndInframeAltering,
                    StartLostAndSpliceRegion  
                    }
impl FromStr for MutationType
{
    type Err=();
    fn from_str(input_str:&str)->Result<MutationType,()>
    {
        match input_str
        {
            "missense"=>Ok(MutationType::MisSense),
            "*missense"=>Ok(MutationType::SMisSense),
            "frameshift"=>Ok(MutationType::FrameShift),
            "*frameshift"=>Ok(MutationType::SFrameShift),
            "inframe_insertion"=>Ok(MutationType::InframeInsertion),
            "*inframe_insertion"=>Ok(MutationType::SInframeInsertion),
            "inframe_deletion"=>Ok(MutationType::InframeDeletion),
            "*inframe_deletion"=>Ok(MutationType::SInframeDeletion),
            "stop_gained"=>Ok(MutationType::StopGained),
            "stop_lost"=>Ok(MutationType::StopLost),
            "*missense&inframe_altering"=>Ok(MutationType::SMisSenseAndInframeAltering),
            "*frameshift&stop_retained"=>Ok(MutationType::SFrameShiftAndStopRetained),
            "*stop_gained&inframe_altering"=>Ok(MutationType::SStopGainedAndInframeAltering),
            "frameshift&stop_retained"=>Ok(MutationType::FrameShiftAndStopRetained),
            "inframe_deletion&stop_retained"=>Ok(MutationType::InframeDeletionAndStopRetained),
            "inframe_insertion&stop_retained"=>Ok(MutationType::InframeInsertionAndStopRetained),
            "stop_gained&inframe_altering"=>Ok(MutationType::StopGainedAndInframeAltering),
            "start_lost"=>Ok(MutationType::StartLost),
            "*stop_gained"=>Ok(MutationType::SStopGained),
            "stop_lost&frameshift"=>Ok(MutationType::StopLostAndFrameShift),
            "missense&inframe_altering"=>Ok(MutationType::MissenseAndInframeAltering),
            "start_lost&splice_region"=>Ok(MutationType::StartLostAndSpliceRegion),
            _=>Err(())
        }
    }
}
/// An enum that store and classify the type of mutated sequence into Sequences which contain string for example KL or NOP, these are mostly
/// associated with missense mutations and infreamce insertions, EndSequences are sequences that ends with * at the end, most commently seen 
/// with frameshift mutations. Lastly, NoSeq is an option used manily to represent sequences that are only composite of *, for example with stop-gained 
/// and stop losts. It derives the Debug, Clone, PartialEq trait and impelement FromStr trait 
/// ## Examples
///``` 
/// use std::str::FromStr;
/// use ppg_rust::data_structures::mutation_ds::MutatedString;
/// let cases=vec!["KLM","NOP*","*",""].iter().map(|case| case.to_string()).collect::<Vec<String>>();
/// assert_eq!(MutatedString::Sequence(cases[0].clone()),MutatedString::from_str(&cases[0]).unwrap());
/// assert_eq!(MutatedString::EndSequence(cases[1].clone()),MutatedString::from_str(&cases[1]).unwrap());
/// assert_eq!(MutatedString::NotSeq,MutatedString::from_str(&cases[2]).unwrap());
///``` 
#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]
pub enum MutatedString
{
    Sequence(String),
    EndSequence(String),
    NotSeq 
}
impl FromStr for MutatedString
{
    type Err=(); 
    fn from_str(input_str:&str)->Result<MutatedString, ()>
    {
        if input_str.is_empty()
        {
            return Err(())
        }
        else if input_str=="*"
        {
            return Ok(MutatedString::NotSeq);
        }
        else if input_str.matches('*').count()!=0 
        {
            return Ok(MutatedString::EndSequence(input_str.to_string()));
        }
        else
        {
            return Ok(MutatedString::Sequence(input_str.to_string()))
        }       
    }
}
/// A struct to store Information related to an amino acid mutation, the four fields stored in the struct are
/// 1. **ref_aa_position** which store the starting position of the mutation in the *reference* sequence 
/// 2. **mut_aa_position** which stores the starting position of the mutation in the *mutated* sequence 
/// 3. **ref_aa** a *MutatedString* instance storing the reference amino acid sequence at the mutational site 
/// 4. **mut_aa** a *MutatedString* instance storing the mutated amino acid sequence at the mutational site
#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]
pub struct MutationInfo
{
    pub ref_aa_position:u16,
    pub mut_aa_position:u16,
    pub ref_aa:MutatedString,
    pub mut_aa:MutatedString,
}

impl MutationInfo
{
    /// A function to create a new MutationInfo instance 
    /// ## Parameters 
    /// 1. ref_aa_position an int, representing the starting position of the mutation in the *reference* sequence 
    /// 2. mut_aa_position an int, representing the starting position of the mutation in the *mutated* sequence 
    /// 3. ref_aa a *MutatedString*, representing the reference amino acid sequence at the mutational site 
    /// 4. mut_aa a *MutatedString*, representing the mutated amino acid sequence at the mutational site
    /// ## Examples 
    ///``` 
    /// use ppg_rust::data_structures::mutation_ds::MutationInfo; 
    /// let ref_pos=32;
    /// let mut_pos=32;
    /// let ref_seq="*".to_string();
    /// let mut_seq="KLM*".to_string();
    /// let eg_case= MutationInfo::new(ref_pos,mut_pos,ref_seq,mut_seq); 
    /// println!("The example has the following structure {:#?}",eg_case); // uses pretty print, notice the numbers are 0-indexed
    ///``` 
    pub fn new(ref_aa_position:u16, mut_aa_position:u16,ref_aa:String,mut_aa:String)->MutationInfo
    {
        MutationInfo
        {
            ref_aa_position:ref_aa_position-1, // rest the index to be 0-indexed
            mut_aa_position:mut_aa_position-1,
            ref_aa:MutatedString::from_str(&ref_aa).unwrap(),
            mut_aa:MutatedString::from_str(&mut_aa).unwrap(),
        }
    }
}
/// An abstract representation for a mutation that is composite mainly of 4 componants 
/// 1. transcrit_name a *String* containing the transcript name 
/// 2. len an i16 int containg the  length of the mutation
/// 3. mut_type  a *MutationType* enum coding for the mutational type the mutational type s
/// 4. mut_info a *MutationInfo* struct summarizing all the mutational info 
///``` 
///``` 
use serde::{Deserialize, Serialize};
#[derive(Debug,Clone,Eq,Serialize,Deserialize)]
pub struct Mutation
{
    pub transcrit_name:String,
    pub mut_type:MutationType,
    pub mut_info:MutationInfo
}
impl Mutation
{
    pub fn new(info_vec:Result<Vec<String>,String>)->Result<Mutation,String>
    {
        let info_vec=match info_vec
        {
            Ok(res)=>res,
            Err(err_msg)=>return Err(format!("Failed to parse the mutation with the following error: {}",err_msg))
        };
        if info_vec.len()!=3
        {
            return Err(format!("Info_vec must be of size 3, however, your input is of size {}",info_vec.len()));
        }
        let mut_type=match  MutationType::from_str(&info_vec[0])
        {
            Ok(mut_type)=>mut_type,
            Err(_)=>
            {
                return Err(format!("The provided mutation: {} is not supported",&info_vec[0]));
            }
        };
        let mut_info= match text_parser::parse_amino_acid_field(&info_vec[2])
        {
            Ok(info_field)=>info_field,
            Err(err_msg)=>
            {
                return Err(format!("Parsing the provided info field: {} failed with the following error message : {}", &info_vec[2], err_msg));
            }
        };
        Ok(Mutation{mut_type:mut_type,mut_info:mut_info,transcrit_name:info_vec[1].clone()})
    }
}
impl Ord for Mutation
{
    fn cmp(&self, other:&Self)->Ordering
    {
        self.mut_info.ref_aa_position.cmp(&other.mut_info.ref_aa_position)
    }
}
impl PartialOrd for Mutation
{
    fn partial_cmp(&self, other:&Self)->Option<Ordering>
    {
        Some(self.cmp(other))
    }
}
impl PartialEq for Mutation
{
    fn eq(&self, other:&Self)->bool
    {
        self.mut_info.mut_aa_position==other.mut_info.mut_aa_position
    }
}
// unit testing the mutational module 
///``` 
///```
#[cfg(test)]
mod test_mutationds
{
    use super::*;
    #[test]
    fn test_mutation_type_enum()
    {
        let test_cases=vec!["missense","*missense","frameshift","*frameshift",
        "inframe_insertion","*inframe_insertion","inframe_deletion","*inframe_deletion","stop_gained", "stop_lost"];
        assert_eq!(MutationType::MisSense, MutationType::from_str(&test_cases[0]).unwrap());
        assert_eq!(MutationType::SMisSense, MutationType::from_str(&test_cases[1]).unwrap());
        assert_eq!(MutationType::FrameShift, MutationType::from_str(&test_cases[2]).unwrap());
        assert_eq!(MutationType::SFrameShift, MutationType::from_str(&test_cases[3]).unwrap());
        assert_eq!(MutationType::InframeInsertion, MutationType::from_str(&test_cases[4]).unwrap());
        assert_eq!(MutationType::SInframeInsertion, MutationType::from_str(&test_cases[5]).unwrap());
        assert_eq!(MutationType::InframeDeletion, MutationType::from_str(&test_cases[6]).unwrap());
        assert_eq!(MutationType::SInframeDeletion, MutationType::from_str(&test_cases[7]).unwrap());
        assert_eq!(MutationType::StopGained, MutationType::from_str(&test_cases[8]).unwrap());
        assert_eq!(MutationType::StopLost, MutationType::from_str(&test_cases[9]).unwrap());
    }
    #[test]
    fn test_mutated_string()
    {
        let test_cases=vec!["KLM","NOP*","*",""].iter().map(|case| case.to_string()).collect::<Vec<String>>();
        assert_eq!(MutatedString::Sequence(test_cases[0].clone()),MutatedString::from_str(&test_cases[0]).unwrap());
        assert_eq!(MutatedString::EndSequence(test_cases[1].clone()),MutatedString::from_str(&test_cases[1]).unwrap());
        assert_eq!(MutatedString::NotSeq,MutatedString::from_str(&test_cases[2]).unwrap());
        match  MutatedString::from_str(&test_cases[3])
        {
        
            Ok(_)=>panic!("Code should failed"),
            _=>()
        }
    }
    #[test]
    fn test_mutation_constructor()
    {
        // define a test-case
        let test_case=vec!["stop_gained".to_string(),"ENST00000484547".to_string(), "32Q>32*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        // assert that it produce the correct results 
        assert_eq!(MutationType::StopGained,test_mutation.mut_type);
        assert_eq!("ENST00000484547".to_string(),test_mutation.transcrit_name);
        assert_eq!(MutationInfo::new(32, 32, "Q".to_string(), "*".to_string()),test_mutation.mut_info);
    } 
    #[test]
    fn test_mutation_bad_input1()->Result<(),String>
    {
        // define a test-case
        let test_case=vec!["stop_gained".to_string(),"ENST00000484547".to_string(), "32Q>32*".to_string(), "Gene".to_string()];
        let test_mutation=Mutation::new(test_case);
        // check that the function failed 
        match test_mutation
        {
            Ok(mutation)=>Err(format!("Test should have failed, however, it results a mutation {:#?}",mutation)),
            Err(_)=>Ok(())
        }    
    }
    #[test]
    fn test_mutation_bad_input2()->Result<(),String>
    {
        // define a test-case
        let test_case=vec!["stop_gainedd".to_string(),"ENST00000484547".to_string(), "32Q>32*".to_string()];
        let test_mutation=Mutation::new(test_case);
        // check that the function failed 
        match test_mutation
        {
            Ok(mutation)=>Err(format!("Test should have failed, however, it results a mutation {:#?}",mutation)),
            Err(_)=>Ok(())
        }    
    }
    #[test]
    fn test_mutation_bad_input3()->Result<(),String>
    {
        // define a test-case
        let test_case=vec!["stop_gainedd".to_string(),"ENST00000484547".to_string(), "32Q32*".to_string()];
        let test_mutation=Mutation::new(test_case);
        // check that the function failed 
        match test_mutation
        {
            Ok(mutation)=>Err(format!("Test should have failed, however, it results a mutation {:#?}",mutation)),
            Err(_)=>Ok(())
        }    
    }
    #[test]
    fn test_mutation_bad_input4()->Result<(),String>
    {
        // define a test-case
        let test_case=vec!["stop_gainedd".to_string(),"ENST00000484547".to_string(), "".to_string()];
        let test_mutation=Mutation::new(test_case);
        // check that the function failed 
        match test_mutation
        {
            Ok(mutation)=>Err(format!("Test should have failed, however, it results a mutation {:#?}",mutation)),
            Err(_)=>Ok(())
        }    
    }
}



