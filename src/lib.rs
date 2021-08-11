/// # Project Description 
///  The crate contains all the function and modules that were utilized to build PPGG (https://github.com/ikmb/ppg)
/// The crate is composite of 5 main modules:
/// 1. Readers which provide a collection of function for reading Fasta and VCF files 
/// 2. Data_structures which is the major engine of the crate, the different data structures provides 
/// a Wide array of struct to abstract and simplify the analysis of genetic data
/// 3. Parts provides a high-level constructs that are build ontop of other parts of the library 
/// 4. Binders provides support for binding the code with C and CUDA codes through FFI 
pub mod readers;
pub mod data_structures;
pub mod parts; 
pub mod functions;
pub mod writers; 
pub mod binders; 


