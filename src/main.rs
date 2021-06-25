use ppgg_rust::{data_structures::Map::IntMap, parts::{cli,io,exec}, writers};
use std::path::{self, Path}; 
use ppgg_rust::writers::write_intmap2json; 
use chrono::{DateTime, Utc};
fn main()
{
    let args = cli::ParsedInput::new(cli::parse_command_line());
    if args.is_verbose
    {
        println!("Reading and loading the VCF file, starting time is: {}",Utc::now())
    }
    let vec_int_repr=io::parse_vcf(Path::new(&args.path2vcf)).unwrap();
    if args.is_verbose
    {
        println!("VCF file have been parsed and encoded into a vector of intermediate representations, finished at: {}",Utc::now()); 
        println!("Loading the Reference file, starting time is: {}",Utc::now()); 
    }
    let ref_seq=io::read_fasta(Path::new(&args.path2fasta)); 
    if args.write_i_map
    {
        println!("Writing the intermediate representation map, starting at: {}", Utc::now());
        let mut base_dir=args.res_path.clone();
        base_dir.push_str("/int_maps");
        let write_path=Path::new(&base_dir); 
        match std::fs::create_dir(write_path)
        {
            Err(err_msg)=>panic!("Creating a directory to write the intermediate representation failed, with the following error: {}",err_msg),
            _=>()
        }; 
        write_intmap2json(write_path,&vec_int_repr).unwrap(); 
    }
    if args.compute_state
    {
        println!("Computing and writing the stats, starting at: {}", Utc::now()); 
        io::compute_and_write_summary(Path::new(&args.res_path), &vec_int_repr); 
        println!("Computing and writing the stats, finished at: {}", Utc::now()); 
        print!("Generating personalized genomes: starting at: {}", Utc::now());
    }
    let vec_per_genomes= exec::execute(vec_int_repr, args.engine.clone(), &ref_seq);
    if args.is_verbose
    {
        println!("Personalized proteomes have been generated, finished at: {}", Utc::now());
    }
    if args.is_verbose
    {
        println!("Write the generated results, starting at: {}", Utc::now())
    }
    io::write_personalized_genomes(vec_per_genomes, args.engine, args.path2vcf);
    if args.is_verbose
    {
        println!("Execution finished at: {}", Utc::now());
    } 
}