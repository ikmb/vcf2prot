use clap::{Arg, App, ArgMatches};
use core::panic;
use std::{path::Path, str::FromStr};
use crate::data_structures::InternalRep::engines::Engine; 

/// ## Summary 
/// A generic representation for the parsed input parameters used by the executable 
#[derive(Debug,Clone)]
pub struct ParsedInput
{
    pub path2vcf:String, 
    pub path2fasta:String,
    pub res_path:String,
    pub engine:Engine, 
    pub compute_state:bool,
    pub is_verbose:bool,
    pub write_i_map:bool,
    pub write_all:bool,
    pub write_compressed:bool,
    pub write_single_thread:bool 
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
            Some(engine)=> 
            {
                let engine=Engine::from_str(engine).unwrap(); 
                match engine
                {
                    Engine::MT | Engine::ST =>engine,
                    Engine::GPU=> panic!("The current version is a CPU-only version with a single-thread (st) and multi-thread (mt) versions only,\
                     however, you asked for a GPU engine, which is not supported in this version. check the project web-page at: https://github.com/ikmb/ppg for more details.")
                }
            
            },
            None=>panic!("The value of engine has not been provided")          
        };
        /* write_e_map:bool, write_i_map:bool */
        let compute_state=args.is_present("stats"); 
        let is_verbose=args.is_present("verbose");
        let write_i_map=args.is_present("write_int_map"); 
        let write_all=args.is_present("write_all_proteins"); 
        let write_compressed = args.is_present("write_compressed");
        let write_single_thread = args.is_present("write_single_thread");
        ParsedInput{path2vcf,path2fasta,res_path,engine,compute_state,is_verbose,write_i_map,write_all,write_compressed,write_single_thread}
    }
}

pub fn parse_command_line()->ArgMatches
{
    App::new("Vcf2prot")
    .version("0.1.4")
    .author("Hesham ElAbd <h.elabd@ikmb.uni-kiel.de>")
    .about("A rust binary that takes as input a FASTA file containing the reference proteome and\
     a VCF file containing the consequence calling and apply the mutations of each patient\
     to the reference file to generate a FASTA file per sample containing the personalized proteome of that individual.\
     The behavior of the tool can controlled using environmental variables currently\n\
     The main environmental variable is the 'NO_TEST' variable which is used to disable other tests and allow the code to run without\
     the runtime cost of the QC. Otherwise all tests will be run. You can selectively run some tests by first exporting the environmental\
     flag 'RUN_SELECTED_TEST' and the export any of the following test FLAGS:
     1. DEBUG_GPU => Inspect the input arrays to the GPU are inspected for indexing error.
     2. DEBUG_CPU_EXEC => Inspect the vector of tasks provided to the input CPU execution engine.
     3. DEBUG_TXP=Transcript_ID => This flag exports a transcript id that will be used for debugging.
     4. INSPECT_TXP => If set, after each transcript is translated into instruction an inspection function will be called to check the correctness of translation
     5. INSPECT_INS_GEN => Inspect the translation process from mutations to instructions and print detailed error messages incase inspection failed.
     6. PANIC_INSPECT_ERR => If set the code will panic if inspecting the translation from mutation to instruction failed.
     
     For more details, see the project webpage at: https://github.com/ikmb/ppg")
    .arg(Arg::new("vcf_file")
        .short('f')
        .long("vcf_file")
        .value_name("FILE")
        .about("A VCF File containing the consequences calling for each sample.")
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
    .arg(Arg::new("write_int_map")
        .short('i')
        .long("write_int_map")
        .required(false)
        .takes_value(false)
        .about("Write an intermediate map containing the observed mutation per transcript per patient to sub directory in the provided output\
        directory, the directory has a predefined name of 'int_maps'. Inside the directory a JSON file containing the\
        intermediate map of each patient is written."))      
    .arg(Arg::new("write_all_proteins")
        .short('a')
        .long("write_all_proteins")
        .required(false)
        .takes_value(false)
        .about("An optional control flag to control the writing behavior of Vcf2prot, if set Vcf2prot will write the altered and the non-altered, i.e.\
        reference sequences, to the fasta file of each proband. This might increase the size of the generated files considerably.\
        By default this option is switched off."))
    .arg(Arg::new("write_compressed")
        .short('c')
        .long("write_compressed")
        .required(false)
        .takes_value(false)
        .about("An optional control flag to control the writing behavior of Vcf2prot, if set Vcf2prot will write the generated fasta files as g-zipped\
        files, i.e. with the extension .fasta.gz, this can be used to decrease the disk space needed by the generated files, especially, \
        when generating 1000s of files.By default this option is switched off. "))    
    .arg(Arg::new("write_single_thread")
        .short('w')
        .long("write_single_thread")
        .required(false)
        .takes_value(false)
        .about("An optional control flag to control the writing behavior of Vcf2prot, if set only one thread is used to write all generated fasta files,\
        by default, this is the case with a single thread engine, i.e. g st, however, this parameter can be used to overwrite this parameter and \
        to enable a single threaded writing of files when a multi-threaded or a GPU engines have been used for parsing and generating the sequences. "))       
    .get_matches()
}

pub fn check_test_state()
{
    println!("Checking the test state ...");
    match std::env::var("NO_TEST") // disable all tests, nno  test is enabled, all other variables are turned off 
    {
        Ok(_)=>
        {
            // disable DEBUG_GPU even if it was exported 
            match std::env::var("DEBUG_GPU")
            {
                Ok(_)=>
                {
                    println!("WARRING:: Conflicting environmental variables found. NO_TEST is True while DEBUG_GPU IS ALSO True ...");
                    println!("INFO:: NO_TEST has a higher precedence ... removing DEBUG_GPU from the list of environment variables"); 
                    std::env::remove_var("DEBUG_GPU");
                },
                Err(_)=>()
            }
            // disabling DEBUG_CPU_EXEC
            match std::env::var("DEBUG_CPU_EXEC")
            {
                Ok(_)=>
                {
                    println!("WARRING:: Conflicting environmental variables found. NO_TEST is True while DEBUG_CPU_EXEC IS ALSO True ...");
                    println!("INFO:: NO_TEST has a higher precedence ... removing DEBUG_CPU_EXEC from the list of environment variables"); 
                    std::env::remove_var("DEBUG_CPU_EXEC");
                },
                Err(_)=>()
            }
            // disabling DEBUG_TXP
            match std::env::var("INSPECT_TXP")
            {
                Ok(_)=>
                {
                    println!("WARRING:: Conflicting environmental variables found. NO_TEST is True while INSPECT_TXP IS ALSO True ...");
                    println!("INFO:: NO_TEST has a higher precedence ... removing INSPECT_TXP from the list of environment variables"); 
                    std::env::remove_var("INSPECT_TXP");
                },
                Err(_)=>()
            }
            // disabling INSPECT_INS_GEN
            match std::env::var("INSPECT_INS_GEN")
            {
                Ok(_)=>
                {
                    println!("WARRING:: Conflicting environmental variables found. NO_TEST is True while INSPECT_INS_GEN IS ALSO True ...");
                    println!("INFO:: NO_TEST has a higher precedence ... removing INSPECT_INS_GEN from the list of environment variables"); 
                    std::env::remove_var("INSPECT_INS_GEN");
                },
                Err(_)=>()
            }
            // disabling PANIC_INSPECT_ERR
            match std::env::var("PANIC_INSPECT_ERR")
            {
                Ok(_)=>
                {
                    println!("WARRING:: Conflicting environmental variables found. NO_TEST is True while PANIC_INSPECT_ERR IS ALSO True ...");
                    println!("INFO:: NO_TEST has a higher precedence ... removing PANIC_INSPECT_ERR from the list of environment variables"); 
                    std::env::remove_var("PANIC_INSPECT_ERR");
                },
                Err(_)=>()
            }
        },
        Err(_)=>
        {
            match std::env::var("RUN_SELECTED_TEST")
            {
                Ok(_)=>state_env_var(),
                Err(_)=>
                {
                    println!("Enabling All Quality Control Tests ... ");
                    match std::env::var("DEBUG_GPU")
                    {
                        Ok(_)=>(),
                        Err(_)=>std::env::set_var("DEBUG_GPU","TRUE")
                    };
                    match std::env::var("DEBUG_CPU_EXEC")
                    {
                        Ok(_)=>(),
                        Err(_)=>std::env::set_var("DEBUG_CPU_EXEC","TRUE")
                    };
                    match std::env::var("INSPECT_TXP")
                    {
                        Ok(_)=>(),
                        Err(_)=>std::env::set_var("INSPECT_TXP","TRUE")
                    };
                    match std::env::var("INSPECT_INS_GEN")
                    {
                        Ok(_)=>(),
                        Err(_)=>std::env::set_var("INSPECT_INS_GEN","TRUE")
                    };
                    match std::env::var("PANIC_INSPECT_ERR")
                    {
                        Ok(_)=>(),
                        Err(_)=>std::env::set_var("PANIC_INSPECT_ERR","TRUE")
                    };
                }
            }   
        }
    }
}

pub fn state_env_var()
{
    println!(" State of the environmental variables is: "); 
    match std::env::var("DEBUG_GPU")
    {
        Ok(_)=>println!("DEBUG_GPU ==> is set "),
        Err(_)=>()
    };
    match std::env::var("DEBUG_CPU_EXEC")
    {
        Ok(_)=>println!("DEBUG_CPU_EXEC ==> is set "),
        Err(_)=>()
    };
    match std::env::var("DEBUG_TXP")
    {
        Ok(transcript_id)=>println!("DEBUG_TXP ==> is set to {}",transcript_id),
        Err(_)=>()
    };
    match std::env::var("INSPECT_TXP")
    {
        Ok(_)=>println!("INSPECT_TXP ==> is set"),
        Err(_)=>()
    };
    match std::env::var("INSPECT_INS_GEN")
    {
        Ok(_)=>println!("INSPECT_INS_GEN ==> is set"),
        Err(_)=>()
    };
    match std::env::var("PANIC_INSPECT_ERR")
    {
        Ok(_)=>println!("PANIC_INSPECT_ERR ==> is set"),
        Err(_)=>()
    };
}