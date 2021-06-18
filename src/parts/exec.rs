use std::collections::HashMap;

use crate::data_structures::InternalRep::engines::Engine; 
use crate::data_structures::Map::IntMap; 
use crate::data_structures::InternalRep::proband_instructions::ProbandInstruction; 
use crate::data_structures::InternalRep::personalized_genome::PersonalizedGenome; 
use rayon::prelude::*; 

pub fn execute(vec_int_repr:Vec<IntMap>, exec_engine:Engine, ref_seq:&HashMap<String,String>)->Vec<PersonalizedGenome>
{
    
    match exec_engine
    {
        Engine::ST=>
        {
            vec_int_repr.into_iter()
            .map(|proband_map|ProbandInstruction::from_intmap(proband_map, exec_engine.clone(),ref_seq))
            .map(|probandMap|PersonalizedGenome::from_proband_instruction(probandMap,exec_engine.clone(),ref_seq))
            .collect::<Vec<PersonalizedGenome>>()
        },
        Engine::MT=>
        {
            vec_int_repr.into_par_iter()
            .map(|proband_map|ProbandInstruction::from_intmap(proband_map, exec_engine.clone(),ref_seq))
            .map(|probandMap|PersonalizedGenome::from_proband_instruction(probandMap,exec_engine.clone(),ref_seq))
            .collect::<Vec<PersonalizedGenome>>()
        },
        _=>panic!("Unknown engine: {:#?}",exec_engine)
    }
}