use ppgg::parts::{cli,io,exec};
use std::path::{Path, PathBuf}; 
use ppgg::writers::write_intmap2json; 
use chrono::Utc;
fn main()
{
    let args = cli::ParsedInput::new(cli::parse_command_line());

    cli::check_test_state(); // print the state of environmental variables 
    if args.is_verbose
    {
        println!("Reading and loading the VCF file, starting time is: {}",Utc::now())
    }
    let vec_int_repr=io::parse_vcf(Path::new(&args.path2vcf),args.engine.clone()).unwrap();
    if args.is_verbose
    {
        println!("VCF file have been parsed and encoded into a vector of intermediate representations, finished at: {}",Utc::now()); 
        println!("Loading the Reference file, starting time is: {}",Utc::now()); 
    }
    let ref_seq=io::read_fasta(Path::new(&args.path2fasta),args.engine.clone()); 
    if args.write_i_map
    {
        println!("Writing the intermediate representation map, starting at: {}", Utc::now());
        let mut pathbuf=PathBuf::from(&args.res_path.clone());
        pathbuf.push("int_maps"); 
        let write_path=Path::new(&pathbuf); 
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
        println!("Generating personalized genomes: starting at: {}", Utc::now());
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
    io::write_personalized_genomes(vec_per_genomes, args.engine, args.res_path,
         args.write_single_thread.clone(),args.write_all.clone(),
         args.write_compressed.clone(), &ref_seq);
    if args.is_verbose
    {
        println!("Execution finished at: {}", Utc::now());
    } 
}