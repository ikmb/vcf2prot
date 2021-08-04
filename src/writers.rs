use std::path::{Path, PathBuf}; 
use std::collections::HashMap; 
use crate::data_structures::Constants;
use crate::data_structures::Map;
use serde_json; 
use std::io::Write;
use std::fs::{File,create_dir};
/// ## Summary 
/// Write the provided earlymap representation into a json file, the function create a directory and write 
/// a JSON file per patient in the directory, the function returns an error if the directory already exists.
pub fn write_earlymap2json(path2write:&Path, vec_earlymap: &Vec<Map::EarlyMap> )->Result<(),String>
{
    match create_dir(path2write)
    {
        Ok(_)=>(),
        Err(err_msg)=>return Err(format!("Creating the output directory failed because of: {} ",err_msg))
    };
    for e_map in vec_earlymap.iter()
    {
        let mut temp_path=path2write.clone().to_path_buf(); 
        temp_path.push(e_map.get_proband_name());
        temp_path.set_extension("json"); 
        let writer= match File::create(temp_path.as_path())
        {
            Ok(file)=>file,
            Err(err_msg)=>return Err(format!("Creating a file for {} failed with the following error message",err_msg))
        };
        serde_json::to_writer(writer, e_map).unwrap();
    }
    Ok(())
}
/// ## Summary 
/// Write the provided intermediate representation into a json file, the function create a directory and write 
/// a JSON file per patient in the directory, the function returns an error if the directory already exists 
/// ## Example 
/// ´´´
/// let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
/// write_intmap2json(Path::new("test_data/test_writer"),&int_map_test).unwrap();
///```
pub fn write_intmap2json(path2write:&Path, vec_intmap: &Vec<Map::IntMap> )->Result<(),String>
{
    match create_dir(path2write)
    {
        Ok(_)=>(),
        Err(err_msg)=>return Err(format!("Creating the output directory failed because of: {} ",err_msg))
    };
    for i_map in vec_intmap.iter()
    {
        let mut temp_path=path2write.clone().to_path_buf(); 
        temp_path.push(i_map.get_name());
        temp_path.set_extension("json"); 
        let writer= match File::create(temp_path.as_path())
        {
            Ok(file)=>file,
            Err(err_msg)=>return Err(format!("Creating a file for {} failed with the following error message",err_msg))
        };
        serde_json::to_writer(writer, i_map).unwrap(); 
    }
    Ok(())
}
/// ## Summary 
/// Write the generated number of mutations per proband to a file 
/// ##Example 
///```rust 
/// let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
/// let test_case=summary::compute_number_mutation_per_proband(&int_map_test); 
/// write_num_number_mutation_per_proband(&Path::new("test_data/number_mut_per_proband.tsv"), test_case).unwrap();
///```
pub fn write_num_number_mutation_per_proband(path2file:&Path,stats_table:HashMap<String,u64>)->Result<(),String>
{
    // set the path 2 buffer 
    let mut pathbuf=PathBuf::from(path2file); 
    pathbuf.push("number_of_mutations_per_proband"); 
    pathbuf.set_extension("tsv"); 
    // create the file
    let mut file_handle=match File::create(pathbuf) 
    {
        Ok(file)=>file,
        Err(err_msg)=>return Err(format!("Creating the file: {:#?} failed due to the following error: {}",path2file, err_msg))    
    };
    // write the file 
    write!(&mut file_handle,"Proband Name \t Number of mutations\n").unwrap();
    for (key,state) in stats_table.iter()
    {
        write!(&mut file_handle,"{},\t{}\n", key, state).unwrap(); 
    }
    Ok(())
}
/// write a TSV table containing the number of mutation per probands 
/// ## Example
///```rust
/// let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
/// let test_case=summary::compute_number_mutation_per_proband(&int_map_test); 
/// write_num_number_mutation_per_proband(&Path::new("test_data/number_mut_per_proband.tsv"), test_case).unwrap();
///```
pub fn write_type_mutations_per_patient(path2file:&Path,stats_table:HashMap<String,Vec<u64>>)->Result<(),String>
{
    // create the path to store the files 
    let mut pathbuf=PathBuf::from(path2file); 
    pathbuf.push("type_of_mutations_per_patient"); 
    pathbuf.set_extension("tsv");
    // create a mutable handle to write the results 
    let mut file_handle= match File::create(&pathbuf) 
    {
        Ok(file)=>file,
        Err(err_msg)=>return Err(format!("Creating the file: {:#?} failed due to the following error: {}", pathbuf, err_msg))    
    };
    // write the files 
    write!(&mut file_handle,"Proband Name\t").unwrap();
    for mutation in Constants::SUP_TYPE.iter()
    {
        write!(&mut file_handle,"{}\t",mutation).unwrap();
    }
    for (key,stats) in stats_table.iter()
    {
        write!(&mut file_handle, "{}\t",key).unwrap();
        for stat in stats
        {
            write!(&mut file_handle,"{}\t",stat).unwrap(); 
        } 
    }
    Ok(())
}
/// write a TSV table containing number of mutations per transcript 
/// ## Example 
///```rust
/// let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
/// let test_case=summary::compute_number_of_mutations_per_transcript(&int_map_test); 
/// write_number_of_mutations_per_transcript(&Path::new("test_data/type_mutation_per_proband.tsv"), test_case).unwrap();
///```
pub fn write_number_of_mutations_per_transcript(path2file:&Path,stats_table:HashMap<String,u64>)->Result<(),String>
{
    // Create the path to store the files 
    let mut pathbuf=PathBuf::from(path2file); 
    pathbuf.push("number_of_mutations_per_transcript"); 
    pathbuf.set_extension("tsv");
    // create a file handle
    let mut file_handle= match File::create(&pathbuf) 
    {
        Ok(file)=>file,
        Err(err_msg)=>return Err(format!("Creating the file: {:#?} failed due to the following error: {}",pathbuf, err_msg))    
    };
    write!(&mut file_handle,"Transcript Name \t Number of mutations\n").unwrap();
    for (key,state) in stats_table.iter()
    {
        write!(&mut file_handle,"{},\t{}\n", key, state).unwrap(); 
    }
    Ok(())
}

#[cfg(test)]
pub mod test_json_parsing
{
    use super::*; 
    use crate::parts::io::parse_vcf; 
    use crate::functions::summary;  
    #[test]
    fn test_intmap2json()
    {
        let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
        write_intmap2json(Path::new("test_data/test_writer"),&int_map_test).unwrap();
    }
    #[test]
    fn test_num_number_mutation_per_proband()
    {
        let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
        let test_case=summary::compute_number_mutation_per_proband(&int_map_test); 
        write_num_number_mutation_per_proband(&Path::new("test_data/number_mut_per_proband.tsv"), test_case).unwrap();
    }
    #[test]
    fn test_type_mutations_per_patient()
    {
        let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
        let test_case=summary::compute_type_mutations_per_patient(&int_map_test); 
        write_type_mutations_per_patient(&Path::new("test_data/type_mutation_per_proband.tsv"), test_case).unwrap();
    }
    #[test]
    fn test_num_mut_per_transcript()
    {
        let int_map_test=parse_vcf(&Path::new("/Users/heshamelabd/projects/test_data/test_case_int1.vcf")).unwrap();
        let test_case=summary::compute_number_of_mutations_per_transcript(&int_map_test); 
        write_number_of_mutations_per_transcript(&Path::new("test_data/num_mutation_per_transcript.tsv"), test_case).unwrap();
    }
}