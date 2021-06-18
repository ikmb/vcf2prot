use clap::{Arg, App, ArgMatches};
use core::panic;
use std::{path::Path, str::FromStr};
use crate::data_structures::InternalRep::engines::Engine; 

#[derive(Debug)]
pub struct ParsedInput
{
    pub path2vcf:String, 
    pub path2fasta:String,
    pub res_path:String,
    pub engine:Engine, 
    pub compute_state:bool,
    pub is_verbose:bool,
    pub write_e_map:bool,
    pub write_i_map:bool
}
impl ParsedInput
{
    pub fn new(args:ArgMatches)->Self
    {
        // parse the path 2 VCF files‚
        let path2vcf= match args.value_of("vcf_file")
        {
            Some(path2file)=>path2file.to_string(),
            None=>panic!("Path to the VCF file has not been provided")
        }; 
        if !(Path::new(&path2vcf).exists())
        {
            panic!("The provided path to the VCF file: {} does not exists",path2vcf)
        }
        // parse the path 2 fasta file 
        let path2fasta= match args.value_of("fasta_ref")
        {
            Some(path2file)=>path2file.to_string(),
            None=>panic!("Path to the fasta file has not been provided")
        }; 
        if !(Path::new(&path2fasta).exists())
        {
            panic!("The provided path to the fasta file: {} does not exists",path2fasta)
        }
        // check the output directory exists
        let res_path= match args.value_of("output_path")
        {
            Some(path2file)=>path2file.to_string(),
            None=>panic!("Path to the fasta file has not been provided")
        }; 
        if !(Path::new(&res_path).exists())
        {
            panic!("The provided path to write the results: {} does not exists",path2fasta)
        }
        // now store the value of the flags
        let engine= match args.value_of("engine") 
        {
            Some(engine)=> Engine::from_str(engine).unwrap(),
            None=>panic!("The value of engine has not been provided")          
        };
        /* write_e_map:bool, write_i_map:bool */
        let compute_state=args.is_present("stats"); 
        let is_verbose=args.is_present("verbose");
        let write_e_map=args.is_present("write_early_map"); 
        let write_i_map=args.is_present("write_int_map");          
        ParsedInput{path2vcf,path2fasta,res_path,engine,compute_state,is_verbose,write_e_map,write_i_map}
    }
}

pub fn parse_command_line()->ArgMatches
{
    App::new("PPGG")
    .version("0.1")
    .author("Hesham ElAbd <h.elabd@ikmb.uni-kiel.de>")
    .about("A rust binary that takes as input a fasta file containing the reference proteome and\
     a vcf file containing the consequence calling and apply the mutations of each patient\
     to the reference file to generate a fasta file per proband containing the personalized proteome of that individual.")
    .arg(Arg::new("vcf_file")
        .short('f')
        .long("vcf_file")
        .value_name("FILE")
        .about("A VCF File containing the consequences calling for each proband.")
        .required(true))
    .arg(Arg::new("fasta_ref")
        .short('r')
        .long("fasta_ref")
        .value_name("FILE")
        .about("A VCF File containing the reference proteome with transcript id as identifiers and protein sequences as the body.")
        .required(true))
    .arg(Arg::new("output_path")
        .short('o')
        .long("output_path")
        .value_name("PATH")
        .about("The path to a directory where fasta files will be written.")
        .required(true))
    .arg(Arg::new("engine")
        .short('g')
        .long("engine")
        .value_name("VALUE")
        .about("The Execution engine, can be any of three values, 'st' for single thread, 'mt' for multiple threads and 'gpu' for\
         for using GPU accelerators.")
        .required(true))
    .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .required(false)
            .takes_value(false)
            .about("If set, print a verbose output about the program state."))
    .arg(Arg::new("stats")
            .short('s')
            .long("stats")
            .required(false)
            .takes_value(false)
            .about("If set, stats are computed and are written to the output directory along with the fasta file"))
    .arg(Arg::new("write_early_map")
            .short('e')
            .long("write_early_map")
            .required(false)
            .takes_value(false)
            .about("Write an early map containing the observed mutation per patient to sub directory in the provided output\
             directory, the directory has a predefined name of 'early_maps'. Inside the directory a JSON file containing the\
             early map of each patient is written"))
    .arg(Arg::new("write_int_map")
        .short('i')
        .long("write_int_map")
        .required(false)
        .takes_value(false)
        .about("Write an intermediate map containing the observed mutation per transcript per patient to sub directory in the provided output\
        directory, the directory has a predefined name of 'int_maps'. Inside the directory a JSON file containing the\
        intermediate map of each patient is written."))       
    .get_matches()
}

