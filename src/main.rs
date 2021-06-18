use ppgg_rust::parts::{cli,io,exec};
use std::path::Path; 
fn main()
{
    let args = cli::ParsedInput::new(cli::parse_command_line());
    let vec_int_repr=io::parse_vcf(Path::new(&args.path2vcf)).unwrap();
    let ref_seq=io::read_fasta(Path::new(&args.path2fasta)); 
    let vec_per_genomes= exec::execute(vec_int_repr, args.engine.clone(), &ref_seq); 
    io::write_personalized_genomes(vec_per_genomes, args.engine, args.path2vcf);
}