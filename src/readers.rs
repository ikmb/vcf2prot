use std::path::Path; 
use std::fs; 
use rayon::prelude::*;
use std::collections::HashMap; 
use crate::data_structures::{vcf_ds,FastaFile,Constants}; 
/// Building a VCF reader that reads an input VCF file and returns a results enums, 
/// the Ok branch contains the probands name and the VCF records that contain the supported mutations
/// while the Err branch contain an error string message.
///  ## Example 
///``` 
/// use std::path::Path;
/// use ppgg_rust::readers;
/// let path=Path::new("/Users/heshamelabd/projects/test_data/dev_case_long_and_short.vcf");
/// let (probands, records)= match readers::read_vcf(path)
/// {
///    Ok(res)=>res,
///    Err(err_msg)=> panic!("Should not have failed!!, ".to_string())
/// }; 
///``` 
pub fn read_vcf(path2load:&Path)->Result<(vcf_ds::Probands,vcf_ds::VCFRecords),String>
{
    // Read the file
    let mut lines= match vcf_helpers::read_file(path2load)
    {
        Ok(lines)=>lines,
        Err(err_msg)=>return Err(err_msg)
    };
    // Get the proband names  
    let proband_names = match vcf_helpers::get_probands_names(&mut lines)
    {
        Ok(lines)=>lines, 
        Err(err_msg)=>return Err(err_msg)
    };
    // Remove the header file
    lines.retain(|line| !line.starts_with('#')); // remove all lines starting 
    // parse the records for QC
    let records= match vcf_helpers::get_records(lines)
    {
        Ok(records)=>records,
        Err(err_msg)=>return Err(err_msg)
    };
    // return the results 
    Ok((vcf_ds::Probands::new(proband_names),vcf_ds::VCFRecords::new(records)))
}
/// Takes as an input the path to a fasta file and return a FastaFile or an error message 
///  ## Example 
///``` 
/// use ppgg_rust::data_structures::FastaFile; 
/// use ppgg_rust::readers::read_fasta_file; 
/// use std::path::Path; 
/// let path2file=Path::new("test_data/test_fasta_data1.fasta");
/// let fasta_file=read_fasta_file(path2file).unwrap(); 
/// assert_eq!(fasta_file.get_records().len(),3); 
/// assert!(fasta_file.is_in_records(&"seq1".to_string()));
///``` 
pub fn read_fasta_file(path2load:&Path)->Result<FastaFile::FastaFile,String>
{
    let lines=match vcf_helpers::read_file(path2load)
    {
        Ok(res)=>res,
        Err(err_msg)=>return Err(err_msg)
    }; 
    let mut records=HashMap::new(); 
    let mut header=String::with_capacity(100); 
    let mut sequence=String::with_capacity(5000); 
    for line in lines
    {
        if line.starts_with('>')
        {
            let line=line.strip_prefix('>').unwrap();
            if header.is_empty() 
            {
                header.push_str(&line);
            }
            else
            {
                records.insert(header.clone(), sequence.clone()); 
                header.clear();
                sequence.clear(); 
                header.push_str(&line); 
            }
        }
        else
        {
            sequence.push_str(&line);
        }   
    }
    // add the final records
    records.insert(header.clone(), sequence.clone()); 
    // check the records are not empty 
    if records.len()==0
    {
        return Err(String::from("The provided, file does not have valid sequence records, parsing it returned 0 record")); 
    }
    Ok(FastaFile::FastaFile::new(records))
}

pub mod vcf_helpers
{
    use super::*;
    /// The function takes the path to a VCF file as an input and returns a vector of strings as an output, 
    /// each string is a line in the input file. 
    ///### Error
    /// incase reading the file failed, the function returns a String containing the error message
    ///## Example 
    ///``` 
    /// use std::path::Path; 
    /// let path = Path::new("/Users/heshamelabd/projects/test_data/dev_file.vcf"); 
    /// let lines= ppgg_rust::readers::vcf_helpers::read_file(&path).unwrap();
    /// for line in lines
    /// {
    ///     println!("{}",line)
    /// }
    ///``` 
    pub fn read_file(path2load:&Path)->Result<Vec<String>, String>
    {
        let file_string = match fs::read_to_string(&path2load) 
        {
            Ok(file_string)=> file_string,
            Err(err_msg)=>
            {
                return Err(format!("\n Function: readers::vcf_helpers::read_file --> could not read the provided file, the following error\
                 was generatied while reading it:\n {} \n", err_msg));
            }
        };
        if file_string.is_empty()
        {
            return Err("\n Function: readers::vcf_helpers::read_file, the provided file is empty \n".to_string()); 
        }
        Ok(file_string.lines().map(|line| line.to_owned()).collect::<Vec<String>>())
    }
    /// Extract the probands name from the VCF file, return a vector of string contain the probands names
    /// ## Example 
    ///``` 
    /// use std::path::Path; 
    /// use ppgg_rust::readers::vcf_helpers; 
    /// let file_path= Path::new("/Users/heshamelabd/projects/test_data/test_file2.vcf");
    /// let mut results = vcf_helpers::read_file(&file_path).unwrap(); 
    /// let res_vec=vcf_helpers::get_probands_names(&mut results).unwrap(); 
    /// for proband in res_vec.iter() {println!("{}",proband)} // print the probands name 
    ///``` 
    pub fn get_probands_names(lines: &mut Vec<String>)->Result<Vec<String>,String>
    {
        // extract the patient line  
        let results_line = match lines.iter_mut().find(|line | line.starts_with("#CHROM"))
        {
            Some(line)=>
            {
                line // a mutable copy to the line  
            },
            None=>
            {
                return Err("Could not find a header line".to_string());
            }
        };
        if results_line.ends_with('\n') || results_line.ends_with('\t')
        {
            results_line.pop();
        }
        // split the lines using the \t seperator 
        let mut res=results_line.split('\t').map(|field| field.to_string()).collect::<Vec<String>>();
        if res.len() <8
        {
            return Err(format!("The provided file does not contain the minimum number of columns, expected a minimum of 8, found: {}",res.len()));
        }
        // remove the first 7 elements 
        let res_clean=res.drain(9..).collect::<Vec<String>>(); 
        // check that there is at least one patient in the file 
        if res_clean.len()==0
        {
            return Err("The file does not contain any patients!!, after removing the madatory ".to_string());
        }
        Ok(res_clean)
    }
    // the wraper for the parallization using massage passing 
    pub fn get_records(lines:Vec<String>)->Result<Vec<String>,String>
    {
        let  res=lines.par_iter()
                            .filter( |&line| return_if_supported(line))
                            .map( |line| line.to_owned())
                            .collect::<Vec<String>>(); 
        if res.len()==0
        {
            return Err("Could not extract any records from the provided file!!".to_string());
        }
        Ok(res)
    }

    /// A helper function that inspect the input record, i.e. one line in the body of the VCF line\
    /// and return the true if it contains at least one supported mutation, otherwise  it returns false which signals\
    /// that the line should not be skipped 
    /// ## Example 
    ///``` 
    /// use ppgg_rust::readers::vcf_helpers; 
    /// let test_line="7\t193407\t7_193407_C_A\tC\tA\t1495\tPASS\tAF=2.5e-05;AQ=1495;ExcessHet=3.0103;QD=14.24;AC=1;AN=32920;BCSQ=missense|FAM20C|ENST00000313766|protein_coding|+|13I>13F|193236A>T\t0|1:0.461538:6,7:13:PASS:99:209,0,186:.:.:10410".to_string();
    /// assert_eq!(true, vcf_helpers::return_if_supported(&test_line))
    ///``` 
    pub fn return_if_supported(line:&String)->bool
    {
        let info_field=line.split('\t').collect::<Vec<&str>>()[7]; 
        let mut BCSQ_field=info_field.split(';').collect::<Vec<&str>>(); 
        BCSQ_field.retain(|&sub_str|sub_str.starts_with("BCSQ="));
        if BCSQ_field.len()==0
        {
            return false; // BCSQ not there 
        }
        let BCSQ_field_str=BCSQ_field[0].split('=').collect::<Vec<&str>>()[1];
        if BCSQ_field_str.contains(',')
        {
            let possible_effects=BCSQ_field_str.split(',').collect::<Vec<&str>>();
            for effect in possible_effects.iter()
            {
                if is_supported_csq(effect)
                {
                    return true;
                }
            }
            return false;
        }
        else
        {
            if is_supported_csq(BCSQ_field_str)
            {
                return true;
            }
            return false;
        }
    }

    /// A helper function that inspect the input string and return True if it can be interpreted by the program or False otherwise 
    /// ## Example 
    ///``` 
    /// use ppgg_rust::readers::vcf_helpers; 
    /// let test_case1="missense|FAM20C|ENST00000313766|protein_coding|+|13I>13F|193236A>T";
    /// let test_case2="splice_region|FAM20C|ENST00000313766|protein_coding";
    /// let test_case3="@246435";
    /// assert_eq!(vcf_helpers::is_supported_csq(test_case1),true);
    /// assert_eq!(vcf_helpers::is_supported_csq(test_case2),false);
    /// assert_eq!(vcf_helpers::is_supported_csq(test_case3),false);
    ///``` 
    pub fn is_supported_csq(csq_str:&str)->bool
    {
        // check the correct number of sperators 
        let num_match=csq_str.matches('|').count();
        if num_match !=6 
        {
            return false; 
        }
        let mut_type=csq_str.split('|').collect::<Vec<&str>>()[0];
        if Constants::SUP_TYPE.contains(&mut_type)
        {
            return true;
        }
        false
    }
}

#[cfg(test)]
pub mod test_vcf_helpers
{
    use super::{Path, read_vcf, read_fasta_file}; 
    use super::vcf_helpers; 
    #[test]
    fn test_read_file1()->Result<(),String>
    {
        let path = Path::new("/Users/heshamelabd/projects/test_data/dev_file.vcf"); 
        match vcf_helpers::read_file(&path)
        {
            Ok(res)=>
            {   
                //for line in res {println!("{}",line)}
                assert_eq!(res.len(),500);
                Ok(())
            },
            Err(err_msg)=>Err(err_msg)
        }
    }
    #[test]
    fn test_read_file2()->Result<(),()>
    {
        let not_define_path=Path::new("not_present_file.vcf");
        match vcf_helpers::read_file(&not_define_path)
        {
            Ok(_)=>
            {   
                Err(())
            },
            Err(_)=>Ok(())
        }

    }
    #[test]
    fn test_read_file3()->Result<(),()>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/wrong_file.vcf");
        match vcf_helpers::read_file(&path)
        {
            Ok(_)=>Err(()),
            Err(_)=>Ok(())
        }
    }
    #[test]
    fn test_get_proband_names()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_file1.vcf");
        let mut results = match vcf_helpers::read_file(&path)
        {
            Ok(res)=>
            {
                res
            },
            Err(err_msg)=>return Err(err_msg)
        }; 
        match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Ok(())
        }; 
        Err("Function did not algin with the expected behaviour of returning an error, because the file is empty".to_string())
    }
    #[test]
    fn test_get_proband_names2()->Result<(),()>
    {
        let case_cor_res="KIEL_ADC00143_0219294502".to_string();
        let file_path= Path::new("/Users/heshamelabd/projects/test_data/test_file2.vcf");
        let mut results = match vcf_helpers::read_file(&file_path)
        {
            Ok(res)=>
            {
                res
            },
            Err(_)=>return Err(())
        }; 
        let res_vec=match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Err(())
        }; 
        assert_eq!(res_vec.len(),1);
        assert_eq!(res_vec[0],case_cor_res);
        Ok(())
    }
    #[test]
    fn test_get_proband_names3()->Result<(),()>
    {
        let res=16460; // number of patient 
        let file_path= Path::new("/Users/heshamelabd/projects/test_data/test_case3.vcf");
        let mut results = match vcf_helpers::read_file(&file_path)
        {
            Ok(res)=>
            {
                res
            },
            Err(_)=>return Err(())
        }; 
        let res_vec=match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Err(())
        }; 
        assert_eq!(res_vec.len(),res); 
        Ok(())
    }
    #[test]
    fn test_get_proband_names4()->Result<(),()>
    {
        let res=16460; // number of patient 
        let file_path= Path::new("/Users/heshamelabd/projects/test_data/test_case4.vcf");
        let mut results = match vcf_helpers::read_file(&file_path)
        {
            Ok(res)=>
            {
                res
            },
            Err(_)=>return Err(())
        }; 
        let res_vec=match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Err(())
        }; 
        assert_eq!(res_vec.len(),res); 
        Ok(())
    }
    #[test]
    fn test_get_proband_names5()->Result<(),()>
    {
        let res=16460; // number of patient 
        let file_path= Path::new("/Users/heshamelabd/projects/test_data/dev_file.vcf");
        let mut results = match vcf_helpers::read_file(&file_path)
        {
            Ok(res)=>
            {
                res
            },
            Err(_)=>return Err(())
        }; 
        let res_vec=match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Err(())
        }; 
        assert_eq!(res_vec.len(),res); 
        Ok(())
    }
    #[test]
    fn test_get_proband_names6()->Result<(),()>
    {
        let first_line="##fileformat=VCFv4.2".to_string();
        let file_path= Path::new("/Users/heshamelabd/projects/test_data/dev_file.vcf");
        let mut results = match vcf_helpers::read_file(&file_path)
        {
            Ok(res)=>
            {
                res
            },
            Err(_)=>return Err(())
        }; 
        match vcf_helpers::get_probands_names(&mut results)
        {
            Ok(res)=>res,
            Err(_)=>return Err(())
        }; 
        assert_eq!(results[0],first_line); 
        assert_eq!(results.len(),500);
        if ! results[499].starts_with('7')
        {
            return Err(())
        }
        Ok(())
    }
    #[test]
    fn test_is_supported()
    {
        let test_case1="intron|FAM20C||protein_coding";
        let test_case2="synonymous|FAM20C|ENST00000313766|protein_coding|+|435Y|256705C>T";
        let test_case3="missense|FAM20C|ENST00000313766|protein_coding|+|440R>440C|256718C>T";
        let test_case4="non_coding|FAM20C||retained_intron";
        let test_case5="splice_region|FAM20C|ENST00000313766|protein_coding";
        let test_case6="@246435";
        let test_case7="missense|FAM20C|ENST00000313766|protein_coding|+|13I>13F|193236A>T";
        let test_case8="*missense|FAM20C|ENST00000313766|protein_coding|+|564N>564D|259915A>G";
        let test_case9="inframe_insertion|FOXL3|ENST00000506382|protein_coding|+|125Y>125YR|291161CCGGCGG>CCGGCGGCGG";
        let test_case10="*inframe_insertion|TMEM184A|ENST00000297477|protein_coding|-|393S>393GS|1547017A>AGCC";
        let test_case11="inframe_deletion|FAM20C|ENST00000313766|protein_coding|+|72AASS>72A|193413GCCGCCTCCT>G";
        let test_case12="*inframe_deletion|FAM20C|ENST00000313766|protein_coding|+|72AASS>72A|193413GCCGCCTCCT>G";
        let test_case13="frameshift|FOXL3|ENST00000506382|protein_coding|+|226VGLHFWTM*>226VDSTFGQC|291463TG>T";
        let test_case14="*frameshift|FOXL3|ENST00000506382|protein_coding|+|226VGLHFWTM*>226VDSTFGQC|291463TG>T";
        let test_case15="stop_gained|FAM20C|ENST00000313766|protein_coding|+|115E>115*|193542G>T";
        let test_case16="stop_gained|FAM20C|ENST00000313766|protein_coding|+|115E>115*|193542G>T";
        let test_case17="stop_lost|FOXL3|ENST00000506382|protein_coding|+|232TM*>232T|291482AATGTG>A"; 
        // check that the function returned the correct results 
        assert_eq!(vcf_helpers::is_supported_csq(test_case1),false);
        assert_eq!(vcf_helpers::is_supported_csq(test_case2),false);
        assert_eq!(vcf_helpers::is_supported_csq(test_case3),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case4),false);
        assert_eq!(vcf_helpers::is_supported_csq(test_case5),false);
        assert_eq!(vcf_helpers::is_supported_csq(test_case6),false);
        assert_eq!(vcf_helpers::is_supported_csq(test_case7),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case8),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case9),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case10),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case11),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case12),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case13),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case14),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case15),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case16),true);
        assert_eq!(vcf_helpers::is_supported_csq(test_case17),true);
    }
    #[test]
    fn test_return_if_supported_case_1()
    {
        let line="7\t193236\t7_193236_A_T\tA\tT\t1789\tPASS\tAF=2.5e-05;AQ=1789;ExcessHet=3.0103;QD=15.29;AC=1;AN=32920;BCSQ=missense|FAM20C|ENST00000313766|protein_coding|+|13I>13F|193236A>T\tGT:AB:AD:DP:GQ:PL:RNC:SBPV:BCSQ 0|0:.:79,0:79:99:.:.:.:0\t0|0:.:76,0:76:99:.:.:.:0\t0|0:.:84,0:84:99:.:.:.:0".to_string(); 
        assert_eq!(true, vcf_helpers::return_if_supported(&line));
    }
    #[test]
    fn test_return_if_supported_case_2()
    {
        let lines_okay=vec![
            "7\t193407\t7_193407_C_A\tC\tA\t1495\tPASS\tAF=2.5e-05;AQ=1495;ExcessHet=3.0103;QD=14.24;AC=1;AN=32920;BCSQ=missense|FAM20C|ENST00000313766|protein_coding|+|70P>70T|193407C>A\tGT:AB:AD:DP:GQ:PL:RNC:SBPV:BCSQ 0|0:.:87,0:87:99:.:.:.:0\t0|0:.:66,0:66:99:.:.:.:0\t0|0:.:108,0:108:99:.:.:.:0\t0|0:.:68,0:68:99:.:.:.:0".to_string(), 
            "7\t291463\t7_291463_TG_T;7_291463_T_G\tTG\tT,GG\t2121\tPASS\tAF=2.5e-05,2.5e-05;AQ=2121,1818;QD=10.2;AC=1,0;AN=32920;BCSQ=frameshift|FOXL3|ENST00000506382|protein_coding|+|226VGLHFWTM*>226VDSTFGQC|291463TG>T,frameshift|FOXL3|ENST00000510017|NMD|+|40VGLHFWTM*>40VDSTFGQC|291463TG>T\tGT:AB:AD:DP:GQ:RNC:SBPV:BCSQ\t0|0:.:16,0,0:16:39:.:.,0:0\t0|0:.:69,0,0:69:99:.:.,0:0\t0|0:.:60,0,0:60:98:.:.,0:0\t0|0:.:65,0,0:65:99:.:.,0:0".to_string(),
            "7\t12342095\t7_12342095_AC_A\tAC\tA\t1259\tPASS\tAC=4;AF=0.000177;AN=32918;AQ=1259;ExcessHet=22.7857;QD=3.84;BCSQ=3_prime_utr|VWDE|ENST00000521169|NMD,3_prime_utr|VWDE|ENST00000452576|NMD,@12351642,*frameshift|VWDE|ENST00000275358|protein_coding|-|1411QCKPGWYGPTCSTALCDPVCLNGGSCNKPNTCLCPNGFFGEHCQNAFCHPPCKNGGHCMRNNVCVCREGYTGRRFQKSICDPTCMNGGKCVGPSTCSCPSGWSGKRCNTPICLQKCKNGGECIAPSICHCPSSWEGVRCQIPICNPKCLYGGRCIFPNVCSCRTEYSGVKCEKKIQIRRH*>1411HANLAGMDPPVVQLCVTLSASMVVRVISQILASVQMDSLGNTVRMLSVTLPVRMVATA*|12337185A>C+12340410G>A+12342095AC>A,frameshift|VWDE|ENST00000275358|protein_coding|-|1411QCKPGWYGPTCSTALCDPVCLNGGSCNKPNTCLCPNGFFGEHCQNAFCHPPCKNGGHCMRNNVCVCREGYTGRRFQKSICDPTCMNGGKCVGPSTCSCPSGWSGKRCNTPICLQKCKNGGECIAPSICHCPSSWEGVRCQIPICNPKCLYGGRCIFPNVCSCRTEYSGVKCEKKIQIRRH*>1411HANLAGMDPPVVQLCATLSASMVVRVISQILASVQMDSLGNTVRMLSVTLPVRMVATA*|12337185A>C+12342095AC>A,frameshift|VWDE|ENST00000614403|protein_coding|-|865QCKPGWYGPTCSTALCDPVCLNGGSCNKPNTCLCPNGFFGEHCQNAFCHPPCKNGGHCMRNNVCVCREGYTGRRFQKRHL*>865HANLAGMDPPVVQLCVTLSASMVVRVISQILASVQMDSLGNTVRMLSVTLPVRMVATA*|12337185A>C+12340410G>A+12342095AC>A,frameshift|VWDE|ENST00000614403|protein_coding|-|865QCKPGWYGPTCSTALCDPVCLNGGSCNKPNTCLCPNGFFGEHCQNAFCHPPCKNGGHCMRNNVCVCREGYTGRRFQKRHL*>865HANLAGMDPPVVQLCATLSASMVVRVISQILASVQMDSLGNTVRMLSVTLPVRMVATA*|12337185A>C+12342095AC>A\tGT:AB:AD:DP:FT:GQ:PL:RNC:SBPV:BCSQ\t0|0:.:45,0:45:PASS:82:.:.:.:0\t0|0:.:48,0:48:PASS:85:.:.:.:0\t0|0:.:47,0:47:PASS:84:.:.:.:0\t0|0:.:57,0:57:PASS:95:.:.:.:0".to_string()
        ]; 
        for line in lines_okay
        {
            assert_eq!(true, vcf_helpers::return_if_supported(&line));
        }
        let lines_bad=vec![
            "7\t2250000\t7_2250000_A_C\tA\tC\t2001\tPASS\tAF=2.5e-05;AQ=2001;ExcessHet=3.0103;QD=24.11;AC=1;AN=32920;BCSQ=splice_region|NUDT1|ENST00000356714|protein_coding,splice_region|NUDT1|ENST00000397046|protein_coding,splice_region|NUDT1|ENST00000397048|protein_coding,splice_region|NUDT1|ENST00000343985|protein_coding,splice_region|NUDT1|ENST00000339737|protein_coding,splice_region|NUDT1|ENST00000397049|protein_coding,@2249999\tGT:AB:AD:DP:GQ:PL:RNC:SBPV:BCSQ 0|0:.:15,0:15:37:.:.:.:0\t0|0:.:19,0:19:44:.:.:.:0\t0|0:.:16,0:16:39:.:.:.:0\t0|0:.:33,0:33:67:.:.:.:0".to_string(),
            "7\t4785683 7_4785683_C_T\tC\tT\t3000\tPASS\tAC=2231;AF=0.069;AN=32918;AQ=3000;ExcessHet=160;QD=0.02;BCSQ=splice_region|AP5Z1|ENST00000649063|protein_coding,splice_region|AP5Z1|ENST00000650310|NMD,splice_region|AP5Z1|ENST00000648925|NMD,splice_region|AP5Z1|ENST00000649315|NMD,splice_region&3_prime_utr|AP5Z1|ENST00000647984|NMD,synonymous|AP5Z1|ENST00000649315|NMD|+|113H|4785683C>T,@4785616,*synonymous|AP5Z1|ENST00000650310|NMD|+|377H|4785683C>T,*synonymous|AP5Z1|ENST00000648925|NMD|+|377H|4785683C>T,*synonymous|AP5Z1|ENST00000649063|protein_coding|+|377H|4785683C>T\tGT:AB:AD:DP:FT:GQ:PL:RNC:SBPV:BCSQ\t0|0:.:43,0:43:PASS:80:.:.:.:0\t0|0:.:45,0:45:PASS:82:.:.:.:0\t0|0:.:15,0:15:PASS:37:.:.:.:0\t0|0:.:40,0:40:PASS:76:.:.:.:0".to_string(),
            "7\t22976212\t7_22976212_T_C\tT\tC\t1445\tPASS\tAC=10836;AF=0.339;AN=31922;AQ=1445;ExcessHet=160;QD=0.01;BCSQ=splice_region|FAM126A|ENST00000432176|protein_coding,splice_region|FAM126A|ENST00000440481|protein_coding,splice_region|FAM126A|ENST00000409923|protein_coding,*synonymous|FAM126A|ENST00000409923|protein_coding|-|208S|22976212T>C,@22976776,*synonymous|FAM126A|ENST00000432176|protein_coding|-|208S|22976212T>C,*synonymous|FAM126A|ENST00000440481|protein_coding|-|259S|22976212T>C\tGT:AB:AD:DP:FT:GQ:PL:RNC:SBPV:BCSQ\t0|1:0.555556:5,5:9:PASS:99:127,0,127:.:.:10410\t0|1:0.631579:12,7:19:PASS:99:209,0,324:..:0.6906:10410\t0|0:.:20,0:20:PASS:46:.:.:.:0\t0|1:0.461538:6,7:13:PASS:99:209,0,186:.:.:10410".to_string()
        ];
        for line in lines_bad
        {
            assert_eq!(false, vcf_helpers::return_if_supported(&line));
        }
    }
    #[test]
    fn test_get_records()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case5.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());} 
            
        };
        match vcf_helpers::get_records(records)
        {
            Ok(_)=>Err("The function should have failed".to_string()),
            Err(_)=>Ok(())
        }
    }
    #[test]
    fn test_get_records2()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case6.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),1);
                Ok(())
            },
            Err(_)=>Ok(())
        }
    }
    #[test]
    fn test_get_records3()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case7.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),28);
                Ok(())
            },
            Err(_)=>Err("Test failed".to_string())
        }
    }
    #[test]
    fn test_get_records4()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case9.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),168);
                Ok(())
            },
            Err(_)=>Err("Test failed".to_string())
        }
    }
    #[test]
    fn test_get_records5()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case10.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),168);
                Ok(())
            },
            Err(_)=>Err("Test failed".to_string())
        } // test_case10.vcf
    }
    #[test]
    fn test_get_records6()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case11.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),3);
                Ok(())
            },
            Err(_)=>Err("Test failed".to_string())
        } // test_case10.vcf
    }
    #[test]
    fn test_get_records7()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/test_case12.vcf");
        let records= match vcf_helpers::read_file(path)
        {
            Ok(res)=>res,
            Err(_)=>{return Err("Something went wrong!".to_string());}   
        };
        match vcf_helpers::get_records(records)
        {
            Ok(parsed_records)=>
            {
                assert_eq!(parsed_records.len(),3);
                Ok(())
            },
            Err(_)=>Err("Test failed".to_string())
        } // test_case10.vcf
    }
    #[test]
    fn test_read_vcf()->Result<(),String>
    {
        let path=Path::new("/Users/heshamelabd/projects/test_data/dev_file.vcf");
        match read_vcf(path)
        {
            Ok(_)=>Ok(()),
            Err(_)=>Err("Failed !!1".to_string())
        }
    }
    #[test]
    fn test_read_vcf2()
    {
        let cases=vec!["KIEL_ADC00143_0219294502",	"KIEL_ADC00167_0219294499",
        "KIEL_ADC00144_0219294489","KIEL_ADC00168_0219294492",	"KIEL_ADC00159_0219294500",	"KIEL_ADC00147_0219294549",	"KIEL_ADC00132_0219294536",
            "KIEL_ADC00133_0219294527",	"KIEL_ADC00130_021929456", "KIEL_ADC00134_0219294512",	"KIEL_ADC00033_0219294479"]; 

        let results = cases.iter().map(|line| line.to_string()).collect::<Vec<String>>(); 

        let path=Path::new("/Users/heshamelabd/projects/test_data/dev_case_long_and_short.vcf");
        let (probands, _)= match read_vcf(path)
        {
            Ok(res)=>res,
            Err(_)=> panic!("Failed !!1")
        }; 
        assert_eq!(results.len(),probands.get_num_probands());

    }
    #[test]
    fn test_fasta_reader1()->Result<(),String>
    {
        let test_case=Path::new("test_data/test_fasta_data1.fasta");
        let records = match read_fasta_file(test_case)
        {
            Ok(res)=>res,
            Err(_)=>return Err("Error in reading the file".to_string())
        };
        // check with a results query 
        let res1="testseq1".to_string(); 
        let res2="testseq2".to_string();
        let res3="testseq3".to_string();
        // assert equal the results 
        assert_eq!(records.get_record(&"seq1".to_string())?,&res1);
        assert_eq!(records.get_record(&"seq2".to_string())?,&res2);
        assert_eq!(records.get_record(&"seq3".to_string())?,&res3);
        Ok(())
    }
    #[test]
    fn test_fasta_reader2()->Result<(),String>
    {
        let test_case=Path::new("test_data/test_fasta_data2.fasta");
        let records = match read_fasta_file(test_case)
        {
            Ok(res)=>res,
            Err(_)=>return Err("Error in reading the file".to_string())
        };
        // check with a results query 
        let res1="testseq1testseq1.1testseq1.2".to_string(); 
        let res2="testseq2".to_string();
        let res3="testseq3testseq3.1testseq3.2".to_string();
        // assert equal the results 
        assert_eq!(records.get_record(&"seq1".to_string())?,&res1);
        assert_eq!(records.get_record(&"seq2".to_string())?,&res2);
        assert_eq!(records.get_record(&"seq3".to_string())?,&res3);
        Ok(())
    }
}







