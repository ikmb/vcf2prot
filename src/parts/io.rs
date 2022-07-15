// load the libraries and crates 
use std::collections::HashMap;
use std::path::Path; 
use rayon::prelude::*; 
use crate::data_structures::InternalRep::engines::Engine;
use crate::data_structures::InternalRep::personalized_genome::PersonalizedGenome;
use crate::readers; 
use crate::data_structures::Map::{self, IntMap}; 
use crate::functions::vcf_tools; 
use crate::parts::exec; 
use crate::writers;
/// ## Summary  
/// Parsing a VCF file and return a result object containing a vector of internal representations
pub fn parse_vcf(path2load:&Path, engine:Engine)->Result<Vec<Map::IntMap>,String>
{
    // Get the proband name 
    let (probands,records)=match readers::read_vcf(path2load, engine.clone()) // clone the engine which is a cheap enum so we can use it later 
    {
        Ok(res)=>res,
        Err(err_msg)=>return Err(format!(" reading the file failed: \n {} \n, formatting the string failed",err_msg))
    }; 
    // Get an early map from the generate probands and records 
    let vec_early_map=vcf_tools::get_early_map(probands, records, engine.clone());    
    // generate an intermediate map 
    Ok(vcf_tools::early_to_intermediate_repr(vec_early_map,engine.clone()))
}
/// ## Summary 
/// Read a fasta file and return a hashmap with sequence id as keys and sequences as values 
pub fn read_fasta(path2load:&Path,engine:Engine)->HashMap<String,String>
{
    readers::read_fasta_file(path2load,engine).unwrap().consume_and_get_hash_map()
}
/// ## Summary 
/// Write the personalized genomes as fasta files to the disk 
pub fn write_personalized_genomes(mut vec_genomes:Vec<PersonalizedGenome>, exec_engines:Engine, output_dir:String,
    use_single_thread:bool, write_all:bool, write_compressed:bool, ref_seq:&HashMap<String,String>)
{
    // this parameter has precedence over the engine and it forces the writing to be carried out in a single threaded manner
    if use_single_thread
    {
        vec_genomes.iter()
            .for_each(|genome|genome.write(&output_dir,&write_all,&write_compressed,&ref_seq).unwrap())
    }
    // if the use_single_thread is not there, then we fallback to the engine guided execution
    match exec_engines
    {
        Engine::ST=>
        {
            vec_genomes.iter()
            .for_each(|genome|genome.write(&output_dir,&write_all,&write_compressed,&ref_seq).unwrap())
        },
        Engine::MT | Engine::GPU=>
        {
            vec_genomes.par_iter_mut()
            .for_each(|genome|genome.write(&output_dir,&write_all,&write_compressed,&ref_seq).unwrap())
        }
    }
}
/// ## Summary 
/// A wrapper function for computing and writing the summary results 
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
