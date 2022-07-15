use crate::data_structures::{InternalRep::engines::Engine, Map::{EarlyMap, IntMap}, vcf_ds::{AltTranscript, Probands, VCFRecords}}; 
use crate::functions::text_parser; 
use rayon::prelude::*;


/// ## Summary 
/// create a vector of early maps on parallele using probands name and a mutable a collection 
/// of VCF records 
pub fn get_early_map(probands:Probands, mut records:VCFRecords, engine:Engine)->Vec<EarlyMap>
{
    // the consequence vector per patient 
    let mutation_per_proband=records.get_csq_per_patient(probands.get_num_probands(),engine.clone()); 

    (probands.get_probands(),mutation_per_proband).into_par_iter()
    .map(|(proband,(vec_mut_one,vec_mut_two))|{
        EarlyMap::new(proband,vec_mut_one,vec_mut_two)})
    .collect::<Vec<EarlyMap>>()
}
/// ## Summary 
/// Process a collection of early maps to a collection of Intermediate maps on Parallel.
pub fn early_to_intermediate_repr(mut vec_of_early_maps:Vec<EarlyMap>,engine:Engine)->Vec<IntMap>
{
    match engine 
    {
        Engine::ST=>
        {
            vec_of_early_maps.iter_mut()
                    .map(|early_map| build_int_map_from_early(early_map))
                    .collect::<Vec<IntMap>>()
        },
        Engine::MT | Engine::GPU =>
        {
            vec_of_early_maps.par_iter_mut()
                    .map(|early_map| build_int_map_from_early(early_map))
                    .collect::<Vec<IntMap>>()
        }
    }
}
/// ## Summary 
/// Build an intermediate map instance, IntMap from an early map instance 
pub fn build_int_map_from_early(early_map:&EarlyMap)->IntMap
{
    // get the map of each mutations in the file 
    let (mutations1,mutations2)=early_map.get_mutations_ref();
    // get the map of each mutations in the file 
    let alt_transcripts1=group_muts_per_transcript(mutations1); 
    let alt_transcripts2=group_muts_per_transcript(mutations2); 
    IntMap::new(early_map.get_proband_name().clone(),alt_transcripts1,alt_transcripts2)
}
/// ## Summary 
/// Group all mutations in each transcript to a vector of AltTranscript, where each element in the generated transcript
/// contain all the mutations observed in a single transcript.
/// ## Example
///```
/// use ppgg_rust::functions::vcf_tools::group_muts_per_transcript; 
/// let mutations=vec!["*missense|MAD1L1|Transcript1|protein_coding|-|1R>1H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript1|protein_coding|-|10R>10H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript2|protein_coding|-|100R>100H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript2|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript3|protein_coding|-|18R>18H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript3|protein_coding|-|1993R>1993H|1936821C>T".to_string(),
///        ];
///        let results=group_muts_per_transcript(&mutations);
///        // check correct length 
///        assert_eq!(results.len(),3);
///        // check the correct name 
///        assert_eq!(results[0].name,"Transcript1");
///        assert_eq!(results[1].name,"Transcript2");
///        assert_eq!(results[2].name,"Transcript3");
///        // check the correct mapping 
///        // 1. results: Transcript 1
///        assert_eq!(results[0].get_alts()[0].mut_info.ref_aa_position,0);
///        assert_eq!(results[0].get_alts()[1].mut_info.ref_aa_position,9);
///        // 2. results: Transcripts 2
///        assert_eq!(results[1].get_alts()[0].mut_info.ref_aa_position,99);
///        assert_eq!(results[1].get_alts()[1].mut_info.ref_aa_position,999);
///        // 3. results: Transcripts 3
///        assert_eq!(results[2].get_alts()[0].mut_info.ref_aa_position,17);
///        assert_eq!(results[2].get_alts()[1].mut_info.ref_aa_position,1992);
///
///```
pub fn group_muts_per_transcript(vec_mut:&Vec<String>)->Vec<AltTranscript>
{
    let mut res=Vec::new(); 
    // define the unique transcripts
    //------------------------------
    for transcript in get_unique_transcript(vec_mut)
    {
        let muts_in_transcript=vec_mut.iter()
                                        .filter(|&file|file.contains(&transcript))
                                        .map(|input_string| input_string.clone())
                                        .collect::<Vec<String>>(); 
        res.push(AltTranscript::new(transcript.clone(), muts_in_transcript))
    }
    res
}
/// ## Summary 
/// Extract the set of uniuqe transcripts in a collection of mutations 
/// ##Example 
///```
/// use ppgg_rust::functions::vcf_tools::get_unique_transcript; 
/// let mutations=vec!["*missense|MAD1L1|Transcript1|protein_coding|-|1R>1H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript1|protein_coding|-|10R>10H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript2|protein_coding|-|100R>100H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript2|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript3|protein_coding|-|100R>100H|1936821C>T".to_string(),
///                "*missense|MAD1L1|Transcript3|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
///        ];
///        let results=get_unique_transcript(&mutations);
///        assert_eq!(results.len(),3);
///        assert_eq!(results[0],"Transcript1");
///        assert_eq!(results[1],"Transcript2");
///        assert_eq!(results[2],"Transcript3");
///
///```
pub fn get_unique_transcript(vec_mut:&Vec<String>)->Vec<String>
{
    let mut uniuqe_muts=vec_mut
                                .iter()
                                .filter_map(|field|
                                    {
                                        match text_parser::split_csq_string(field)
                                        {
                                            Ok(res)=>Some(res[1].clone()),
                                            Err(_)=>None,
                                        }
                                    })
                                .collect::<Vec<String>>(); 
    uniuqe_muts.sort(); 
    uniuqe_muts.dedup(); 
    uniuqe_muts

}
 
#[cfg(test)]
mod test_vcf_tools_function
{
    use super::*; 
    #[test]
    pub fn test_get_early_map()
    {
        let probands=vec![
                    "Proband1".to_string(),
                    "Proband2".to_string(),
                    "Proband3".to_string(),
                    "Proband4".to_string()]; 
        let mutations=vec![
            ("mutation1_1,mutation1_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),"mutation1_2,mutation1_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()),
            ("mutation2_1,mutation2_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),"mutation2_2,mutation2_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()),
            ("mutation3_1,mutation3_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),"mutation3_2,mutation3_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()),
            ("mutation4_1,mutation4_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),"mutation4_2,mutation4_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>())
            ]; 
        let results=(probands,mutations).into_par_iter()
        .map(|(proband,(vec_mut_one,vec_mut_two))|{
            EarlyMap::new(proband,vec_mut_one,vec_mut_two)})
        .collect::<Vec<EarlyMap>>(); 
        // write the assert statment 
        assert_eq!(results[0].get_proband_name(),"Proband1");
        assert_eq!(results[1].get_proband_name(),"Proband2");
        assert_eq!(results[2].get_proband_name(),"Proband3");
        assert_eq!(results[3].get_proband_name(),"Proband4");
        // check that the mutations have been inserted 
        assert_eq!(results[0].get_mutations_ref(),(&"mutation1_1,mutation1_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),
                                                &"mutation1_2,mutation1_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()));
        assert_eq!(results[1].get_mutations_ref(),(&"mutation2_1,mutation2_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),
                                                &"mutation2_2,mutation2_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()));
        assert_eq!(results[2].get_mutations_ref(),(&"mutation3_1,mutation3_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),
                                                &"mutation3_2,mutation3_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()));
        assert_eq!(results[3].get_mutations_ref(),(&"mutation4_1,mutation4_3".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>(),
                                                &"mutation4_2,mutation4_4".split(",").map(|elem|elem.to_string()).collect::<Vec<String>>()));
    }
    #[test]
    pub fn test_get_unique_transcript()
    {
        let mutations=vec!["*missense|MAD1L1|Transcript1|protein_coding|-|1R>1H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript1|protein_coding|-|10R>10H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript2|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript2|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript3|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript3|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
        ];
        let results=get_unique_transcript(&mutations);
        assert_eq!(results.len(),3);
        assert_eq!(results[0],"Transcript1");
        assert_eq!(results[1],"Transcript2");
        assert_eq!(results[2],"Transcript3");
    }
    #[test]
    fn test_group_muts_per_transcript()
    {
        let mutations=vec!["*missense|MAD1L1|Transcript1|protein_coding|-|1R>1H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript1|protein_coding|-|10R>10H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript2|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript2|protein_coding|-|1000R>1000H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript3|protein_coding|-|18R>18H|1936821C>T".to_string(),
                "*missense|MAD1L1|Transcript3|protein_coding|-|1993R>1993H|1936821C>T".to_string(),
        ];
        let results=group_muts_per_transcript(&mutations);
        // check correct length 
        assert_eq!(results.len(),3);
        // check the correct name 
        assert_eq!(results[0].name,"Transcript1");
        assert_eq!(results[1].name,"Transcript2");
        assert_eq!(results[2].name,"Transcript3");
        // check the correct mapping 
        // 1. results: Transcript 1
        assert_eq!(results[0].get_alts()[0].mut_info.ref_aa_position,0);
        assert_eq!(results[0].get_alts()[1].mut_info.ref_aa_position,9);
        // 2. results: Transcripts 2
        assert_eq!(results[1].get_alts()[0].mut_info.ref_aa_position,99);
        assert_eq!(results[1].get_alts()[1].mut_info.ref_aa_position,999);
        // 3. results: Transcripts 3
        assert_eq!(results[2].get_alts()[0].mut_info.ref_aa_position,17);
        assert_eq!(results[2].get_alts()[1].mut_info.ref_aa_position,1992);
    }
}










