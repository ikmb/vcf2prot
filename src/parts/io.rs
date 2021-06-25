use core::panic;
use std::collections::HashMap;
use std::fmt::write;
use std::path::Path; 
use rayon::prelude::*; 
use crate::data_structures::InternalRep::engines::Engine;
use crate::data_structures::InternalRep::personalized_genome::PersonalizedGenome;
use crate::readers; 
use crate::data_structures::Map::{self, IntMap}; 
use crate::functions::vcf_tools; 
use crate::parts::exec; 
use crate::writers; 
/// parsing a VCF file and return a result object containing a vector of internal representations
pub fn parse_vcf(path2load:&Path)->Result<Vec<Map::IntMap>,String>
{
    // Get the proband name 
    let (probands,records)=match readers::read_vcf(path2load)
    {
        Ok(res)=>res,
        Err(err_msg)=>return Err(format!(" reading the file failed: \n {} \n, formatting the string failed",err_msg))
    }; 
    // Get an early map from the generate probands and records 
    let vec_early_map=vcf_tools::get_early_map(probands, records); 
    // generate an intermediate map 
    Ok(vcf_tools::early_to_intermediate_repr(vec_early_map))
}
/// read a fasta file and return a hashmap with sequence id as keys and sequences as values 
pub fn read_fasta(path2load:&Path)->HashMap<String,String>
{
    readers::read_fasta_file(path2load).unwrap().consume_and_get_hash_map()
}
/// write the personalized genomes as fasta files to the disk 
pub fn write_personalized_genomes(mut vec_genomes:Vec<PersonalizedGenome>, exec_engines:Engine, output_dir:String)
{
    match exec_engines
    {
        Engine::ST=>
        {
            vec_genomes.iter()
            .for_each(|genome|genome.write(&output_dir).unwrap())
        },
        Engine::MT=>
        {
            vec_genomes.par_iter_mut()
            .for_each(|genome|genome.write(&output_dir).unwrap())
        }
        _=>panic!("Unknown Engine: {:#?}",exec_engines)
    }
}

pub fn compute_and_write_summary(path2write:&Path, vec_maps:&Vec<IntMap>)
{
    // compute the stats, write the files sequentially 
    let computed_stats=exec::compute_states(&vec_maps); 
    let (mut_per_patient,type_mut_per_patient, num_mut_per_transcript )=(
        computed_stats.num_mutation_per_proband,computed_stats.type_mutation_per_proband,
        computed_stats.number_of_mutations_per_transcript); // get the results as three variables 
    // write the results
    writers::write_num_number_mutation_per_proband(path2write, mut_per_patient).unwrap(); 
    writers::write_type_mutations_per_patient(path2write, type_mut_per_patient).unwrap(); 
    writers::write_number_of_mutations_per_transcript(path2write, num_mut_per_transcript).unwrap(); 
}
