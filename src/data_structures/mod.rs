/// ## Summary
/// The module directory contain different modules and sub-modules arranged as follow:
/// 1. mutation_ds ==> contain mutational parsing functions mainly focusing on consequence calling strings for example 
/// inframe_deletion|HOXB3|ENST00000471459|protein_coding|-|89GG>89G|48551141TCCGCCGCCGCCGCCA>TCCGCCGCCGCCA 
/// 2. vcf_ds ==> contains structs and function for handling and working with VCF files 
/// 3. FastaFile ==> a wrapper for reading Fasta files and for parsing them 
/// 4. InternalRep ==> is a module made from mandy structures and submodule and represent the back bone for SIR based representation of sequences
/// 5. Map ==> contains structures for handling the mapping between probands in the VCF files and there corresponding mutation 
/// 6. MaskDecoder ==> contains the class bitmask decoder 
/// 7. Constants ==> contains constant values used throughput the library 
pub mod mutation_ds;
pub mod vcf_ds; 
pub mod FastaFile;
pub mod InternalRep; 
pub mod Map; 
pub mod MaskDecoder;
pub mod Constants; 
