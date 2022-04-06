/// The module contains function to compute statistical summaries 
use std::collections::HashMap; 
use rayon::prelude::*; 
use crate::data_structures::Map; 
use crate::data_structures::mutation_ds::MutationType; 
use std::str::FromStr;

/// ## Summary
/// Compute the mutational load per patient from an input vector of intermediate representation 
pub fn compute_number_mutation_per_proband(vec_intmaps:&Vec<Map::IntMap>)->HashMap<String,u64>
{
    vec_intmaps.par_iter()
    .map(|int_map|stat_helper::number_mutations_per_proband(int_map))
    .collect::<HashMap<String,u64>>()
}
/// ## Summary
/// Compute the number of each mutational type in the patient from an input vector of intermediate representations 
pub fn compute_type_mutations_per_patient(vec_intmaps:&Vec<Map::IntMap>)->HashMap<String,Vec<u64>>
{
    vec_intmaps.par_iter()
    .map(|int_map|stat_helper::compute_mutation_frequency(int_map))
    .collect::<HashMap<String,Vec<u64>>>()
}
/// ## Summary
/// Compute the number of alterations or mutations per transcript from an input vector of intermediate representations 
pub fn compute_number_of_mutations_per_transcript(vec_intmaps:&Vec<Map::IntMap>)->HashMap<String,u64>
{
    let transcript_names=stat_helper::get_uniuqe_transcript(vec_intmaps);
    transcript_names.into_par_iter()
    .map(|transcript_name| stat_helper::count_in_all_individual(transcript_name,&vec_intmaps))
    .collect::<HashMap<String,u64>>()
}

mod stat_helper
{
    use crate::data_structures::Constants;

    use super::*; 
    pub fn number_mutations_per_proband(int_map:&Map::IntMap)->(String,u64)
    {
        let (mut1,mut2)=int_map.get_mutations_ref();
        let num_mut=mut1.len()+mut2.len(); 
        (int_map.get_name().clone(),num_mut as u64)
    }

    pub fn compute_mutation_frequency(int_map:&Map::IntMap)->(String,Vec<u64>)
    {
        let counts=Constants::SUP_TYPE.par_iter()
        .map(|&mut_type|get_count_per_proband(mut_type,int_map))
        .collect::<Vec<u64>>();        
        (int_map.proband_name.clone(),counts)
    }

    fn get_count_per_proband(mut_type:&str, int_map:&Map::IntMap)->u64
    {
        let mut sum=0; 
        let (mut_h1, mut_h2)= int_map.get_mutations_ref();
        let mut_type=MutationType::from_str(mut_type).unwrap(); 
        for alt in mut_h1.iter()
        {
            for mutation in alt.get_alts().iter()
            {
                if mutation.mut_type==mut_type
                {
                    sum+=1
                }
            }
        }
        for alt in mut_h2.iter()
        {
            for mutation in alt.get_alts().iter()
            {
                if mutation.mut_type==mut_type
                {
                    sum+=1
                }
            }
        }
        sum
    }
    pub fn get_uniuqe_transcript(vec_intmaps:&Vec<Map::IntMap>)->Vec<String>
    {
        let mut results=vec_intmaps.par_iter()
        .map(|map|extract_transcript_from_map(map))
        .flatten()
        .collect::<Vec<String>>(); 
        results.sort();
        results.dedup();
        results
    }
    fn extract_transcript_from_map(intmap:&Map::IntMap)->Vec<String>
    {
        let (mut_h1,mut_h2)=intmap.get_mutations_ref();
        let mut res1=mut_h1.par_iter()
        .map(|alt|alt.name.clone())
        .collect::<Vec<String>>(); 
        let mut res2=mut_h2.par_iter()
        .map(|alt|alt.name.clone())
        .collect::<Vec<String>>(); 
        res1.append(&mut res2); 
        res1
    }
   pub fn count_in_all_individual(transcript_name:String,vec_intmaps:&Vec<Map::IntMap>)->(String,u64)
   {
        let results=vec_intmaps.par_iter()
        .map(|intmap|get_count_in_a_proband(&transcript_name,intmap))
        .collect::<Vec<u64>>(); 
    (transcript_name,results.iter().sum())
   }
   fn get_count_in_a_proband(transcript_name:&String, intmap:&Map::IntMap)->u64
   {
        let (mut_h1,mut_h2)=intmap.get_mutations_ref(); 
        let mut sum=0;
        if mut_h1.iter().any(|alt|alt.name==*transcript_name){sum+=1}
        if mut_h2.iter().any(|alt|alt.name==*transcript_name){sum+=1}
        sum
   }
}
#[cfg(test)]
pub mod test_summary_function
{
    use super::*; 
    use crate::{parts::io::parse_vcf, data_structures::InternalRep::engines::Engine}; 
    fn generate_default_internal_representation()->Vec<Map::IntMap>
    {      
        use std::path::Path; 
        match parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf"),Engine::MT)
        {
            Ok(res)=>res,
            Err(err_msg)=>panic!("{}",err_msg)
        }
    }
    #[test]
    fn test_compute_number_mutation_per_proband()
    {
        let num_mut_per_pat=compute_number_mutation_per_proband(&generate_default_internal_representation());
        println!("{:#?}",num_mut_per_pat); 

    }
    #[test]
    fn test_compute_type_mutations_per_patient()
    {
        let type_mutation_per_proband=compute_type_mutations_per_patient(&generate_default_internal_representation());
        println!("{:#?}",type_mutation_per_proband); 
    }
    #[test]
    fn test_number_of_mutations_per_transcript()
    {
        let num_mut_per_transcript=compute_number_of_mutations_per_transcript(&generate_default_internal_representation());
        println!("{:#?}",num_mut_per_transcript); 
    }
}