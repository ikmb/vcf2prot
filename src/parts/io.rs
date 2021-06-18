use core::panic;
use std::collections::HashMap;
use std::path::Path; 
use rayon::prelude::*; 
use crate::data_structures::InternalRep::engines::Engine;
use crate::data_structures::InternalRep::personalized_genome::PersonalizedGenome;
use crate::readers; 
use crate::data_structures::Map; 
use crate::functions::vcf_tools; 
// parsing VCF file 
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

pub fn read_fasta(path2load:&Path)->HashMap<String,String>
{
    readers::read_fasta_file(path2load).unwrap().consume_and_get_hash_map()
}

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