use std::collections::HashMap;
use crate::functions::summary::*; 
use crate::data_structures::InternalRep::engines::Engine; 
use crate::data_structures::Map::IntMap; 
use crate::data_structures::InternalRep::proband_instructions::ProbandInstruction; 
use crate::data_structures::InternalRep::personalized_genome::PersonalizedGenome; 
use rayon::prelude::*; 
use crossbeam::thread; 

// drive the public functions 
//---------------------------
#[derive(Debug,Clone)]
pub struct StatSummary
{
    pub num_mutation_per_proband:HashMap<String,u64>,
    pub type_mutation_per_proband:HashMap<String,Vec<u64>>,
    pub number_of_mutations_per_transcript:HashMap<String,u64>,
}
/// The executioner for computing and generating a personalized proteome per patient 
pub fn execute(vec_int_repr:Vec<IntMap>, exec_engine:Engine, ref_seq:&HashMap<String,String>)->Vec<PersonalizedGenome>
{
    match exec_engine
    {
        Engine::ST=>
        {
            vec_int_repr.into_iter()
            .map(|proband_map|ProbandInstruction::from_intmap(proband_map, exec_engine.clone(),ref_seq))
            .map(|proband_map|PersonalizedGenome::from_proband_instruction(proband_map,exec_engine.clone(),ref_seq))
            .collect::<Vec<PersonalizedGenome>>()
        },
        Engine::MT | Engine::GPU =>
        {
            vec_int_repr.into_par_iter()
            .map(|proband_map|ProbandInstruction::from_intmap(proband_map, exec_engine.clone(),ref_seq))
            .map(|probandMap|PersonalizedGenome::from_proband_instruction(probandMap,exec_engine.clone(),ref_seq))
            .collect::<Vec<PersonalizedGenome>>()
        }
    }
}
/// A function to compute the state from the vec_maps, it launches 3 threads to compute each metric on parallel
pub fn compute_states(vec_maps:&Vec<IntMap>)->StatSummary
{
    thread::scope(|scope|
    {
        // launch the threads 
        let mutation_per_proband=scope.spawn(|_|compute_number_mutation_per_proband(vec_maps)); 
        let type_mutation_per_proband=scope.spawn(|_|compute_type_mutations_per_patient(vec_maps)); 
        let number_mut_per_transcript=scope.spawn(|_|compute_number_of_mutations_per_transcript(vec_maps)); 
        // wait for the results 
        let mut_per_proband=mutation_per_proband.join().unwrap(); 
        let type_mut_per_proband=type_mutation_per_proband.join().unwrap(); 
        let number_mut_per_transcript=number_mut_per_transcript.join().unwrap(); 
        // return the results 
        StatSummary
        {
            num_mutation_per_proband:mut_per_proband,
            type_mutation_per_proband:type_mut_per_proband,
            number_of_mutations_per_transcript:number_mut_per_transcript,
        }   
    }).unwrap()
}