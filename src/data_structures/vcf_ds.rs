use core::panic;

use rayon::prelude::*; 
use crate::functions::text_parser; 
use crate::data_structures::{MaskDecoder::BitMask,
                            mutation_ds::Mutation
                            };

use super::Constants; 
use serde::{Deserialize, Serialize};
/// An abstraction for a collection of VCF Records, the struct owns the provided vector of strings,
/// where each string is a record from the file.
#[derive(Debug,Clone)]
pub struct VCFRecords
{
    records:Vec<String>,
}
impl VCFRecords
{
    /// Create a new VCFRecords from a vector of strings  where each string represent a line in the VCF file 
    /// ## Example
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::VCFRecords; 
    /// let test_case1=vec![
    ///        "1\t1\t1\t1\t1\t1\t1\t1\tField1\tField1.1\tField1.2\tField1.3\tField1.4\tField1.5\tField1.6\tField1.7".to_string(),
    ///        "1\t1\t1\t1\t1\t1\t1\t1\tField2\tField2.1\tField2.2\tField2.3\tField2.4\tField2.5\tField2.6\tField2.7".to_string(),
    ///        ];
    /// let records=VCFRecords::new(test_case1); 
    ///``` 
    pub fn new(records:Vec<String>)->Self
    {
        VCFRecords{records}
    }
    /// Return a reference 
    /// ## Example
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::VCFRecords; 
    /// let test_case1=vec![
    ///        "1\t1\t1\t1\t1\t1\t1\t1\tField1\tField1.1\tField1.2\tField1.3\tField1.4\tField1.5\tField1.6\tField1.7".to_string(),
    ///        "1\t1\t1\t1\t1\t1\t1\t1\tField2\tField2.1\tField2.2\tField2.3\tField2.4\tField2.5\tField2.6\tField2.7".to_string(),
    ///        ];
    /// let records=VCFRecords::new(test_case1.clone()); 
    /// let results=records.get_records(); 
    /// assert_eq!(results,&test_case1); 
    ///```
    pub fn get_records(&self)->&Vec<String>
    {
        &self.records
    }
    /// Return a vector of String containing the consequences at each genomic position 
    /// ## Example
    ///```
    /// use std::path::Path; 
    /// use ppgg_rust::readers; 
    /// use ppgg_rust::data_structures::vcf_ds::VCFRecords; 
    /// let case_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_case.vcf"); 
    /// let res_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_res.txt");
    /// let records=readers::vcf_helpers::read_file(case_path).unwrap(); 
    /// let results=readers::vcf_helpers::read_file(res_path)
    ///                                .unwrap()
    ///                                .iter()
    ///                                .map(|elem| elem.split("=").map(|elem| elem.to_string()).collect::<Vec<String>>().last().unwrap().clone())
    ///                                .collect::<Vec<String>>();
    /// let vcf_records=VCFRecords::new(records);
    /// assert_eq!(results,vcf_records.get_consequences_vector());
    ///```
    pub fn get_consequences_vector(&self)->Vec<String>
    {
        self.records.par_iter()
                    .map(|rec|rec.split("\t").collect::<Vec<&str>>()[7])
                    .map(|rec| rec.split("BCSQ=").collect::<Vec<&str>>()[1].to_string())
                    .collect::<Vec<String>>()
    }
    /// Returns a 2D Matrix that is implemented as a vector of vector of strings, where all the mutations observed in a patient is 
    /// collected into one vector, i.e. the inner vector, the outer vector represent the collection of all mutations.
    /// ## Example
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::VCFRecords; 
    /// let test_case1=vec![
    ///    "1\t1\t1\t1\t1\t1\t1\t1\tField1\tField1.1\tField1.2\tField1.3\tField1.4\tField1.5\tField1.6\tField1.7".to_string(),
    ///    "1\t1\t1\t1\t1\t1\t1\t1\tField2\tField2.1\tField2.2\tField2.3\tField2.4\tField2.5\tField2.6\tField2.7".to_string(),
    ///    "1\t1\t1\t1\t1\t1\t1\t1\tField3\tField3.1\tField3.2\tField3.3\tField3.4\tField3.5\tField3.6\tField3.7".to_string(),
    ///    "1\t1\t1\t1\t1\t1\t1\t1\tField4\tField4.1\tField4.2\tField4.3\tField4.4\tField4.5\tField4.6\tField4.7".to_string(),
    ///    "1\t1\t1\t1\t1\t1\t1\t1\tField5\tField5.1\tField5.2\tField5.3\tField5.4\tField5.5\tField5.6\tField5.7".to_string(),
    ///    ];
    /// let test_result1=vec![
    ///        "Field1", "Field2", "Field3", "Field4", "Field5"].
    ///        iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
    /// let test_result2=vec![
    ///        "Field1.1", "Field2.1", "Field3.1", "Field4.1", "Field5.1"].
    ///        iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
    /// let test_result3=vec![
    ///            "Field1.2", "Field2.2", "Field3.2", "Field4.2", "Field5.2"].
    ///            iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
    /// let test_result4=vec![
    ///                "Field1.3", "Field2.3", "Field3.3", "Field4.3", "Field5.3"].
    ///                iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
    /// let test_result5=vec![
    ///                    "Field1.4", "Field2.4", "Field3.4", "Field4.4", "Field5.4"].
    ///                    iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
    /// let mut vcf_records=VCFRecords::new(test_case1);
    /// assert_eq!(test_result1,vcf_records.get_patient_fields(8)[0]);
    /// assert_eq!(test_result2,vcf_records.get_patient_fields(8)[1]);
    /// assert_eq!(test_result3,vcf_records.get_patient_fields(8)[2]);
    /// assert_eq!(test_result4,vcf_records.get_patient_fields(8)[3]);
    /// assert_eq!(test_result5,vcf_records.get_patient_fields(8)[4]);
    ///```
    pub fn get_patient_fields(&mut self, num_probands:usize)->Vec<Vec<String>>
    {
        let number_probands=&self.records[0].matches('\t').count()-8;
        let mut res:Vec<Vec<String>>=Vec::with_capacity(number_probands);
        for _ in 0..num_probands
        {
            res.push(Vec::with_capacity(self.records.len())); 
        }        
        for record in self.records.iter_mut()
        {
            // parse records  
            let mut fields=record
                            .split("\t")
                            .map(|field|field.to_string())
                            .collect::<Vec<String>>();
            let fields = fields.drain(9..).collect::<Vec<String>>();
            //println!("fields are: {:#?}",fields);
            //println!("Index is: {}, While number of proband is: {}",fields.len(),num_probands); 
            // push to all vectors on parallele 
            for idx in 0..fields.len()
            {
                res[idx].push(fields[idx].clone())
            }
        }
        res
    }

    pub fn get_csq_per_patient(&mut self,num_probands:usize)->Vec<(Vec<String>,Vec<String>)>
    {
        let consequences=self.get_consequences_vector(); 
        let probands_table=self.get_patient_fields(num_probands);
        // we need to get the consequences of each vector 
        probands_table.par_iter()
        .map(|donor|VCFRecords::decode_back(&consequences,donor))
        .collect::<Vec<(Vec<String>,Vec<String>)>>()   
    }

    pub fn decode_back(consequences:&Vec<String>,proband_fields:&Vec<String>)->(Vec<String>,Vec<String>)
    {
        // get index of each conseuqences 
        let mut bitmasks=proband_fields
                            .par_iter()// now only for amoment 
                            .map(|field| text_parser::get_bit_mask(field))
                            .collect::<Vec<String>>();
        // get a vector of tuples at each position 
        let results=(consequences,&mut bitmasks)
                            .into_par_iter()
                            .map(|(csq,bitmask)|VCFRecords::extract_effects(&csq,bitmask))
                            .filter(|(elem1,elem2)|elem1.len()!=0 || elem2.len()!=0)
                            .collect::<Vec<(Vec<String>,Vec<String>)>>(); 
        
        // unroll the mutation into two vectors one for the first haplotype and one for the second 
        let tuple_1_res=results.par_iter()
                                                            .map(|elem|elem.0.clone())
                                                            .flatten()
                                                            .filter(|csq|Constants::SUP_TYPE.contains(&text_parser::get_type(csq)))
                                                            .collect::<Vec<String>>();
        let tuple_2_res=results.par_iter()
                                                            .map(|elem|elem.1.clone())
                                                            .flatten()
                                                            .filter(|csq|Constants::SUP_TYPE.contains(&text_parser::get_type(csq)))
                                                            .collect::<Vec<String>>(); 
        (tuple_1_res,tuple_2_res)
    }
    /// A helper associated function that recieves as an input a CSQ string and a bit mask, results is: Tuple of size two, first element 
    /// is a vector of strings, while number two is a vector of strings => 
    /// each of them are the effects as a specifc site.
    ///```
    /// use std::path::Path; 
    /// use ppgg_rust::readers;
    /// use ppgg_rust::data_structures::vcf_ds::VCFRecords; 
    /// let case_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_case.vcf"); 
    /// let res_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_res.txt");
    /// let records=readers::vcf_helpers::read_file(case_path).unwrap(); 
    /// let resullts=readers::vcf_helpers::read_file(res_path)
    ///                                    .unwrap()
    ///                                    .iter()
    ///                                    .map(|elem| elem.split("=").map(|elem| elem.to_string()).collect::<Vec<String>>().last().unwrap().clone())
    ///                                    .collect::<Vec<String>>();
    //// let vcf_records=VCFRecords::new(records);
    ///```
    pub fn extract_effects(csq:&String, bitmask:&mut String)->(Vec<String>,Vec<String>)
    {
        let splitted_csq=csq.split(",")
                                     .collect::<Vec<&str>>();
        let (haplotype1,haplotype2)=match BitMask::from_string(bitmask).get_indices()
        {
            Some(vec)=>vec,
            None=> return (Vec::new(),Vec::new())
        };
        let index_haplotype_1=haplotype1.into_iter()
            .map(|idx|splitted_csq[idx as usize].to_string())
            .collect::<Vec<String>>();
        let index_haplotype_2=haplotype2.into_iter()
            .map(|idx|splitted_csq[idx as usize].to_string())
            .collect::<Vec<String>>();
        (index_haplotype_1,index_haplotype_2)
    }
}
/// a struct that acts as a wrapper for vector of string containing the name of probands in the VCF file
#[derive(Debug,Clone)]
pub struct Probands
{
    probands:Vec<String>,
}
impl Probands
{
    /// create a new instance using a vector of strings 
    pub fn new(probands:Vec<String>)->Self
    {
        Probands{probands}
    }
    /// return the number of probands in the instance 
    pub fn get_num_probands(&self)->usize
    {
        self.probands.len()
    }
    /// return the vector of strings containing the name of all probands 
    pub fn get_probands(self)->Vec<String>
    {
        self.probands
    }
}
/// AltTranscript => altered transcript, a transcript that is used as an abstraction for 
/// a collection of mutation occuring in a transcript.
#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct AltTranscript
{
    pub name:String,
    pub alts:Vec<Mutation>
}
 impl AltTranscript 
 {
    /// Create a new instance from a transcript name and a vector of Consequences string
    /// ## Example 
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::AltTranscript; 
    /// let name="ENST00000406869".to_string(); 
    /// let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
    ///    ];
    ///    let results=AltTranscript::new(name, mutations); 
    ///    assert_eq!(results.name,"ENST00000406869".to_string());
    ///    assert_eq!(results.get_alts().len(),4);
    ///```    
    pub fn new(name:String,alts:Vec<String>)->Self
    {
        
        let alts=alts.iter()
            .map(|field|Mutation::new(text_parser::split_csq_string(&field).unwrap()).unwrap())
            .collect::<Vec<Mutation>>(); 
        AltTranscript{name,alts}
    }
    /// create a new instance for a transcript name and a vector of mutations that will be filled later 
    pub fn allocate(name:String)->Self
    {
        let alts:Vec<Mutation>=Vec::new(); 
        AltTranscript{name,alts}
    }
    /// create a new instance for a transcript name and a vector of mutations that have an expected number of mutations 
    pub fn with_capacity(name:String, expected_number:usize)->Self
    {
        let alts:Vec<Mutation>=Vec::with_capacity(expected_number); 
        AltTranscript{name,alts}
    }
    /// add an alteration i.e. a genetic mutation, to the current instance of mutations, 
    /// ## Example 
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::AltTranscript; 
    /// let name="ENST00000406869".to_string(); 
    /// let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
    ///    ];
    ///    let mut results=AltTranscript::new(name, mutations); 
    ///    assert_eq!(results.get_alts().len(),4);
    ///    results.add_altes("*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()); 
    ///    assert_eq!(results.get_alts().len(),5);
    ///```    
    pub fn add_altes(&mut self, alt:String)
    {
        self.alts.push(Mutation::new(text_parser::split_csq_string(&alt).unwrap()).unwrap());
    }
    /// return a reference to the instance vector of mutations 
    pub fn get_alts(&self)->&Vec<Mutation>
    {
        &self.alts
    }
    /// add an alteration i.e. a genetic mutation, to the current instance of mutations, 
    /// ## Example 
    ///```
    /// use ppgg_rust::data_structures::vcf_ds::AltTranscript; 
    /// let name="ENST00000406869".to_string(); 
    ///    let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
    ///            "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
    ///    ];
    ///    let mut results=AltTranscript::new(name, mutations); 
    ///    results.add_altes("*missense|MAD1L1|ENST00000406869|protein_coding|-|200L>1000H|1936821C>T".to_string()); 
    ///    results.sort_alterations(); 
    ///    assert_eq!(results.get_alts()[0].mut_info.ref_aa_position,0);
    ///    assert_eq!(results.get_alts()[1].mut_info.ref_aa_position,9);
    ///    assert_eq!(results.get_alts()[2].mut_info.ref_aa_position,99);
    ///    assert_eq!(results.get_alts()[3].mut_info.ref_aa_position,199);
    ///    assert_eq!(results.get_alts()[4].mut_info.ref_aa_position,999);
    ///```    
    pub fn sort_alterations(&mut self)
    {
        self.alts.sort(); 
    }
 }

#[cfg(test)]
mod test_alt_transcript
{
    use super::*; 
    #[test]
    fn test_new()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
        ];
        let results=AltTranscript::new(name, mutations); 
        assert_eq!(results.name,"ENST00000406869".to_string());
        assert_eq!(results.get_alts().len(),4);
    }
    #[test]
    fn test_add_to_alts()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
        ];
        let mut results=AltTranscript::new(name, mutations); 
        assert_eq!(results.get_alts().len(),4);
        results.add_altes("*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()); 
        assert_eq!(results.get_alts().len(),5);
    }
    #[test]
    fn test_sortting_function()
    {
        let name="ENST00000406869".to_string(); 
        let mutations=vec!["*missense|MAD1L1|ENST00000406869|protein_coding|-|1R>1H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|10R>10H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|100R>100H|1936821C>T".to_string(),
                "*missense|MAD1L1|ENST00000406869|protein_coding|-|1000R>1000H|1936821C>T".to_string()
        ];
        let mut results=AltTranscript::new(name, mutations); 
        results.add_altes("*missense|MAD1L1|ENST00000406869|protein_coding|-|200L>1000H|1936821C>T".to_string()); 
        results.sort_alterations(); 
        assert_eq!(results.get_alts()[0].mut_info.ref_aa_position,0);
        assert_eq!(results.get_alts()[1].mut_info.ref_aa_position,9);
        assert_eq!(results.get_alts()[2].mut_info.ref_aa_position,99);
        assert_eq!(results.get_alts()[3].mut_info.ref_aa_position,199);
        assert_eq!(results.get_alts()[4].mut_info.ref_aa_position,999);
    }
    
}
#[cfg(test)]
mod test_vcf
{
    use super::*; 
    use crate::readers; 
    use std::path::Path;
    #[test]
    fn test_get_consequence_file()
    {
        // Prepare input test and results case 
        //------------------------------------
        let case_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_case.vcf"); 
        let res_path=Path::new("/Users/heshamelabd/projects/test_data/test_f1_mod_test_vcf_res.txt");
        let records=readers::vcf_helpers::read_file(case_path).unwrap(); 
        let resullts=readers::vcf_helpers::read_file(res_path)
                                        .unwrap()
                                        .iter()
                                        .map(|elem| elem.split("=").map(|elem| elem.to_string()).collect::<Vec<String>>().last().unwrap().clone())
                                        .collect::<Vec<String>>();
        // Parse the input and check the correctness 
        //------------------------------------------
        let vcf_records=VCFRecords::new(records);
        assert_eq!(resullts,vcf_records.get_consequences_vector());
    }
    #[test]
    fn test_get_patient_fields1()
    {
        let test_case1=vec![
            "1\t1\t1\t1\t1\t1\t1\t1\t1\tField1\tField1.1\tField1.2\tField1.3\tField1.4\tField1.5\tField1.6\tField1.7".to_string(),
            "1\t1\t1\t1\t1\t1\t1\t1\t1\tField2\tField2.1\tField2.2\tField2.3\tField2.4\tField2.5\tField2.6\tField2.7".to_string(),
            "1\t1\t1\t1\t1\t1\t1\t1\t1\tField3\tField3.1\tField3.2\tField3.3\tField3.4\tField3.5\tField3.6\tField3.7".to_string(),
            "1\t1\t1\t1\t1\t1\t1\t1\t1\tField4\tField4.1\tField4.2\tField4.3\tField4.4\tField4.5\tField4.6\tField4.7".to_string(),
            "1\t1\t1\t1\t1\t1\t1\t1\t1\tField5\tField5.1\tField5.2\tField5.3\tField5.4\tField5.5\tField5.6\tField5.7".to_string(),
            ];
        let test_result1=vec![
                "Field1", "Field2", "Field3", "Field4", "Field5"].
                iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
        let test_result2=vec![
                "Field1.1", "Field2.1", "Field3.1", "Field4.1", "Field5.1"].
                iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
        let test_result3=vec![
                    "Field1.2", "Field2.2", "Field3.2", "Field4.2", "Field5.2"].
                    iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
        let test_result4=vec![
                        "Field1.3", "Field2.3", "Field3.3", "Field4.3", "Field5.3"].
                        iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
        let test_result5=vec![
                            "Field1.4", "Field2.4", "Field3.4", "Field4.4", "Field5.4"].
                            iter().map(|elem |elem.to_string()).collect::<Vec<String>>();
        let mut vcf_records=VCFRecords::new(test_case1);
        assert_eq!(test_result1,vcf_records.get_patient_fields(8)[0]);
        assert_eq!(test_result2,vcf_records.get_patient_fields(8)[1]);
        assert_eq!(test_result3,vcf_records.get_patient_fields(8)[2]);
        assert_eq!(test_result4,vcf_records.get_patient_fields(8)[3]);
        assert_eq!(test_result5,vcf_records.get_patient_fields(8)[4]);
    }
    #[test]
    fn test_extract_effects1()
    {
        let mut test_case="effect1,effect2,effect3,effect4,effect5,effect6".to_string(); 
        let mut bit_mask="1".to_string();
        let csq_map=VCFRecords::extract_effects(&mut test_case,&mut bit_mask);
        assert_eq!(csq_map.0[0],"effect1");
    }
    #[test]
    fn test_extract_effects2()
    {
        let mut test_case="effect1,effect2,effect3,effect4,effect5,effect6".to_string(); 
        let mut bit_mask="3".to_string();
        let csq_map=VCFRecords::extract_effects(&mut test_case,&mut bit_mask);
        assert_eq!(csq_map.0[0],"effect1");
        assert_eq!(csq_map.1[0],"effect1");
    }
   
}








