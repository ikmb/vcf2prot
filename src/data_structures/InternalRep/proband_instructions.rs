use std::collections::HashMap;
use super::{engines::Engine, haplotype_instruction::HaplotypeInstruction};
use crate::data_structures::Map::IntMap;

use serde::{Deserialize, Serialize};
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct ProbandInstruction
{
    pub proband_name:String, 
    pub haplotype1_instruction:HaplotypeInstruction, 
    pub haplotype2_instruction:HaplotypeInstruction
}
impl ProbandInstruction
{
    pub fn new(proband_name:String, haplotype1_instruction:HaplotypeInstruction, 
    haplotype2_instruction:HaplotypeInstruction)->Self
    {
        ProbandInstruction{proband_name,haplotype1_instruction,haplotype2_instruction}
    }
    pub fn from_intmap(mut int_map:IntMap, engine:Engine, ref_seq:&HashMap<String,String>)->Self
    {
        let proband_name=int_map.proband_name.clone();
        let (haplo1_vec,haplo2_vec)=int_map.consume_and_get_vecs(); 
        let h1_t_ins= HaplotypeInstruction::from_vec_t_ins(haplo1_vec, engine.clone(),ref_seq); 
        let h2_t_ins= HaplotypeInstruction::from_vec_t_ins(haplo2_vec, engine.clone(),ref_seq);  
        ProbandInstruction::new(proband_name, h1_t_ins, h2_t_ins)
    }
}