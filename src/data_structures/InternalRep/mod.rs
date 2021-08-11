/// The sub-modules contain functions for generating and manipulating SIR representations
/// 1. instruction ==> for the encoding of a mutation into an instruction 
/// 2. transcript_instructions ==> For a collection of instructions observed in a transcript 
/// 3. haplotype_instruction ==> For a collection of mutated transcripts
/// 4. proband_instructions ==> For a collection of two haplotypes representing all the alterations in a transcript 
/// 5. sequence_tape ==> Write the generated sequence into a fasta file 
/// 6. personalized_genome ==> A wrapper for two sequence-tapes used to represent the alteration in a transcript 
/// 7. task ==> a representation for generation a sequence 
/// 8. gir ==> a representation for generating tasks
pub mod instruction; 
pub mod transcript_instructions;
pub mod haplotype_instruction;  
pub mod proband_instructions; 
pub mod sequence_tape; 
pub mod personalized_genome; 
pub mod task; 
pub mod engines;
pub mod gir; 
