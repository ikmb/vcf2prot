use crate::data_structures::mutation_ds::MutationInfo; 
use crate::data_structures::Constants; 
/// The function takes the consequence string and returned a Result enum either containing an Ok or Err type.
/// # Ok
/// a vector of strings that has three elements first the type of the mutation,
/// secnd the transcript id and third the change in the position and the sequence of the mutated amino acids.
/// # Errors 
/// incase the provided string does not have six pipe sperators
///``` 
/// use ppgg_rust::functions::text_parser::split_csq_string;
/// let example_csq_string="stop_gained|RABGEF1|ENST00000484547|NMD|+|32Q>32*|66771993C>T".to_string();
/// match split_csq_string(&example_csq_string)
/// {
///      Ok(fields) =>
///      {
///        for field in fields {println!("{}",field)};
///      }
///     _ =>()
/// }
///```
pub fn split_csq_string(input_string:&String)->Result<Vec<String>,String>
{
    let num_match=input_string.matches('|').count();
    let res = input_string.split('|').map(|elem| elem.into()).collect::<Vec<String>>(); 
    match num_match
    {
        6=>
        {
            if res[3].as_str()!="protein_coding"
            {
                return Err("Skipping this transcript as it is not a protein coding transcript".to_string())   
            }
            let index:Vec<usize>=vec![0,2,5];
            Ok(index.iter().map(|i| res[*i].clone()).collect::<Vec<String>>())
        }, 
        _=>
        {
            match res[0].as_str()
            {
                "start_lost"=>
                {
                    let mut results=Vec::with_capacity(3); 
                    results.push(res[0].clone());
                    results.push(res[2].clone());
                    results.push("1M>1*".to_string()); 
                    Ok(results)
                },
                _=>
                {
                    println!("In correct number of fields, expected 6, received {} and the input string is: {}, skipping this mutation ...",num_match,input_string); 
                    Err(format!("In correct number of fields, expected 6, received {} and the input string is: {}",num_match,input_string))
                }
            }
        }
    }
}
/// The function takes the mutation amino acid field, e.g. "32Q>32*" and returned a Result enum either containing an Ok or Err type.
/// # Ok
/// a MutationInfo struct containg the position of the mutation in the reference and in the mutated amino acids, along with sequence 
/// representation for the mutated and reference seuqence
/// # Errors 
/// incase parsing the provided sequecne failed, a string coding for the error message will be retrained 
///```rust 
/// use ppgg_rust::functions::text_parser;
/// use ppgg_rust::data_structures::mutation_ds::{MutatedString,MutationInfo};
/// let mut_string="32Q>32*".to_string();
/// let res = text_parser::parse_amino_acid_field(&mut_string).expect("Generating the parse_amino_acid failed");
/// let mut_info=MutationInfo
/// {
///     ref_aa_position:31, // zero-based indexing 
///     mut_aa_position:31, // zero-based indexing 
///     ref_aa:MutatedString::Sequence("Q".to_string()),
///     mut_aa:MutatedString::NotSeq,
/// };
/// assert_eq!(mut_info.ref_aa_position,res.ref_aa_position); 
/// assert_eq!(mut_info.mut_aa_position,res.mut_aa_position); 
/// assert_eq!(MutatedString::NotSeq,res.mut_aa); 
/// assert_eq!(MutatedString::Sequence("Q".to_string()),res.ref_aa);
///```
pub fn parse_amino_acid_field(input_string: &String)->Result<MutationInfo,String>
{
    // split the field into two amino acids 
    let parsed_strings=input_string.split('>').collect::<Vec<&str>>();
    if parsed_strings.len()!=2
    {
        return Err(format!("The psrsed string has a length of: {}, expected only two",parsed_strings.len()));
    }
    // get the position and the reference sequence
    let (ref_pos, ref_seq)=match parse_amino_acid_seq_position(&parsed_strings[0])
    {
        Ok((index,sequence))=>(index,sequence), 
        Err(err_msg)=>
        {
            return Err(format!("\n while extracting the sequence and the position of the reference the following error was encounterred {}",err_msg));
        }
    };
    // get the position and the mutated sequence
    let (mut_pos,mut_seq)= match  parse_amino_acid_seq_position(&parsed_strings[1])
    {
        Ok((index,sequence))=>(index,sequence),
        Err(err_msg)=>
        {
            return Err(format!("\n while extracting the sequence and the position of the mutation the following error was encounterred {}",err_msg));
        }
    };
    Ok(MutationInfo::new(ref_pos,mut_pos,ref_seq,mut_seq))
}
/// The function takes an input string composite of an aminoacid position concatinated with a stirng object ,e.g 35KTEST and returns 
/// the amino acid position as u16 int, in this case it 35, and the string containg the mutation, here it is KTEST.
/// ## Ok
/// a tuple containg the amino acid position as an int and the sequence as a stirng,
/// ## Errors 
/// a string contain the cause of faliure
/// # Example
///``` 
/// let test_example="35KTEST";
/// use ppgg_rust::functions::text_parser::parse_amino_acid_seq_position; 
/// match parse_amino_acid_seq_position(test_example)
/// {
///       Ok((pos,seq))=>println!("The position is: {}, while the sequence is: {}",pos,seq),
///       Err(seq)=>()
/// }
///```
pub fn parse_amino_acid_seq_position(input_seq: &str)->Result<(u16,String),String>
{
    if input_seq.matches('-').count() !=0
    {
        return Err(format!("Input string: {} is invalid, it contains a '-' sign which is not valid for indexing amino acid positions, also it is not avalid amino acid",input_seq));
    }
    let input_as_vec=input_seq.chars().collect::<Vec<char>>();// split the input string into a vector of chars, for example, 32Q -> 3,2,Q;
    let nums=['0','1','2','3','4','5','6','7','8','9']; // valid numbers 
    let position = match input_as_vec.iter().clone().filter(|c|  nums.contains(&c)).collect::<String>().parse() // extract the numbers from the stream 
    {
        Ok(num)=>num,
        Err(err_msg)=>
        {
            return Err(format!("Parsing the input sequence {}, failed with the following error message {}",input_seq,err_msg ));
        }
    };
    let mut sequence = input_as_vec.iter().clone().filter(|c| !nums.contains(&c)).collect::<String>();
    if sequence.is_empty()
    {
        sequence="*".to_string(); 
    }
    Ok((position,sequence)) // get a position, sequence tuple 
}  
/// takes an input patient field and extract the bitmask from the last fields of the stirng, e.g. 1|1:3 =>depending on the input
/// it returns either "" an empty string, representing the reference, 3$ if a single int  bit mask or a trimmed versionof the bitmask if 
/// more than one number are provided as input, for example, 1|1:1234,5,0,0,0 => 1234,5 will be the bitmask 
/// # Example
///``` 
/// use ppgg_rust::functions::text_parser;
/// let mut test_case="0|1:0.432432:16,21:37:PASS:99:634,0,417:..:0.1989:10922"; 
/// let mut results=text_parser::get_bit_mask(&test_case.to_string());
/// assert_eq!(results,"10922$"); 
/// test_case="0|1:0.432432:16,21:37:PASS:99:634,0,417:..:0.1989:10922,14,0,0,0";
/// let mut results=text_parser::get_bit_mask(&test_case.to_string());
/// assert_eq!(results,"10922,14"); 
///```
pub fn get_bit_mask(input_string:&String)->String
{
    // check there is at least one semicolon in the patient fields  
    if input_string.matches(":").count()==0
    {
        return Constants::DEF_CONSEQ.to_string();
    }
    // define the bitmask field 
    let bitmask_field=input_string.split(":")
                        .collect::<Vec<&str>>()[input_string.matches(":").count()]// gets the las field of the patient string 
                        .to_string();
    // get the strings 
    if bitmask_field==".".to_string()
    {
        return Constants::DEF_CONSEQ.to_string();
    }
    if bitmask_field.matches(",").count()==0
    {
        return parse_fields(bitmask_field);
    }
    let bitmask_field=remove_leading_zeros(bitmask_field);

    if bitmask_field==Constants::DEF_CONSEQ
    {
        return bitmask_field;
    }
    if bitmask_field.matches(",").count()==0
    {
        return parse_fields(bitmask_field);
    } 
    bitmask_field
}
/// Parse the input string, by trying to cast it as an i32, if this was a Legitimate operation
/// it return the input string with $ appended at the end, for example, 3-->3$ or 0 becomes 0$,
/// if this fail, for whatever reasons, the function returns  Constants::DEF_CONSEQ, currently equal ""
/// # Example
///``` 
/// use ppgg_rust::functions::text_parser::parse_fields; 
/// let test_example="3";
/// let results=parse_fields(test_example.to_string()); 
/// assert_eq!(results,"3$"); 
///```
pub fn parse_fields(mut fields:String)->String
{
    match &mut fields.parse::<i32>()
    {
        Ok(res)=>
        {
            if *res < 0 
            {
                panic!("An invalid bit mask was encountered: {} .  Most likely an outdated version of csq has been used. Check this commit @ Github for more details: https://github.com/samtools/bcftools/commit/1f1e7667ffc1235f31a82e2093f037338acbb4e7",fields);
            }
            fields.push_str("$");
            fields
        }
        Err(_)=> Constants::DEF_CONSEQ.to_string()
    }
}

/// Trim leading zeros from a bitmask string, e.g. 3,5,0->3,5, this is a helper function used to remove the leadng zeros 
/// # Example
///``` 
/// use ppgg_rust::functions::text_parser::remove_leading_zeros; 
/// let test_example="3,5,0";
/// let results=remove_leading_zeros(test_example.to_string()); 
/// assert_eq!(results,"3,5"); 
///```
pub fn remove_leading_zeros(mut fields:String)->String
{
    let mut split_result=fields.split(",")
                .map(|elem| elem.to_string())
                .collect::<Vec<String>>(); 

    while split_result.len()!=0 && split_result.last().unwrap()=="0"
    {
        split_result.remove(split_result.len()-1);
    }
    //fields=fields.trim_end_matches(|char| char=='0' || char ==',').to_string(); 
    if split_result.len()==0
    {
        return Constants::DEF_CONSEQ.to_string();
    }
    else
    {
        if fields.contains('-'){panic!("An invalid bit mask was encountered: {}. Most likely an outdated version of csq has been used. Check this commit @ Github for more details: https://github.com/samtools/bcftools/commit/1f1e7667ffc1235f31a82e2093f037338acbb4e7",fields,);}
    }
    fields=split_result.join(",");
    fields
}
/// a one-liner function for generating the type of mutation from the consequence string. 
/// ## Example 
///``` 
/// use ppgg_rust::functions::text_parser::get_type; 
/// let test_case="*missense|ITPRID1|ENST00000409210|protein_coding|+|717C>717Y|31643796G>A".to_string(); 
/// assert_eq!(*"*missense", *get_type(&test_case));
///```
#[inline]
pub fn get_type(mut_type:&String)->&str
{
    mut_type.split('|').collect::<Vec<&str>>()[0]
}

#[cfg(test)]
mod test_text_parser
{
    use super::*; // make the content of the file public  
    #[test]
    /// The function test that the function produce the correct results when given a default correct string 
    /// This test will be passed if the function only produce three strings from the provided string 
    fn test_split_csq_string_one()->Result<(),String>
    {
        let test_string="stop_gained|RABGEF1|ENST00000484547|NMD|+|32Q>32*|66771993C>T".to_string();
        match split_csq_string(&test_string)
        {
            Ok(test_string)=>
            {
                let res:Vec<String>=vec!["stop_gained","ENST00000484547","32Q>32*"].iter().map(|&x| x.to_string()).collect();
                assert_eq!(res,test_string);
                Ok(())
            },
            Err(error)=>
            {
                Err(format!("Test have failed, the function failed while it should not have failed, error was {}", error))
            }
        }
    }
    #[test]
    fn test_split_csq_string_two()->Result<(),String>
    {
        let test_string_2="5_prime_utr|RABGEF1|ENST00000437078|protein_coding".to_string(); 
        match split_csq_string(&test_string_2)
        {
            Ok(res)=>
            {
                Err(format!("The function should have returned an Err results, however, it returned {}", res[0]))
            },
            Err(error)=>
            {
                Ok(())
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_field_1()
    {
        let mut_string="32Q>32*".to_string();
        let res = parse_amino_acid_field(&mut_string).expect("Generating the parse_amino_acid failed");
        let mut_info=MutationInfo
        {
            ref_aa_position:31,
            mut_aa_position:31,
            ref_aa:MutatedString::Sequence("Q".to_string()),
            mut_aa:MutatedString::NotSeq,
        };
        assert_eq!(mut_info.ref_aa_position,res.ref_aa_position); 
        assert_eq!(mut_info.mut_aa_position,res.mut_aa_position); 
        assert_eq!(MutatedString::NotSeq,res.mut_aa); 
        assert_eq!(MutatedString::Sequence("Q".to_string()),res.ref_aa); 
    }
    #[test]
    fn test_parse_amino_acid_field_2()
    {
        let mut_string="32QK>32*".to_string();
        let res = parse_amino_acid_field(&mut_string).expect("Generating the parse_amino_acid failed");
        let mut_info=MutationInfo
        {
            ref_aa_position:31,
            mut_aa_position:31,
            ref_aa:MutatedString::Sequence("QK".to_string()),
            mut_aa:MutatedString::NotSeq,
        };
        assert_eq!(mut_info.ref_aa_position,res.ref_aa_position); 
        assert_eq!(mut_info.mut_aa_position,res.mut_aa_position); 
        assert_eq!(mut_info.ref_aa,res.ref_aa);
        assert_eq!(mut_info.mut_aa,res.mut_aa); 
    }
    #[test]
    fn test_parse_amino_acid_field_3()
    {
        let mut_string="32QK>32NMKLOPLMNBJK*".to_string();
        let res = parse_amino_acid_field(&mut_string).expect("Generating the parse_amino_acid failed");
        let mut_info=MutationInfo
        {
            ref_aa_position:31,
            mut_aa_position:31,
            ref_aa:MutatedString::Sequence("QK".to_string()),
            mut_aa:MutatedString::EndSequence("NMKLOPLMNBJK*".to_string()),
        };
        assert_eq!(mut_info.ref_aa_position,res.ref_aa_position); 
        assert_eq!(mut_info.mut_aa_position,res.mut_aa_position); 
        assert_eq!(mut_info.ref_aa,res.ref_aa);
        assert_eq!(mut_info.mut_aa,res.mut_aa);
    }
    #[test]
    fn test_parse_amino_acid_field_4()
    {
        let mut_string="32*>32NMKLOPLMNBJK*".to_string();
        let res = parse_amino_acid_field(&mut_string).expect("Generating the parse_amino_acid failed");
        let mut_info=MutationInfo
        {
            ref_aa_position:31,
            mut_aa_position:31,
            ref_aa:MutatedString::NotSeq,
            mut_aa:MutatedString::EndSequence("NMKLOPLMNBJK*".to_string()),
        };
        assert_eq!(mut_info.ref_aa_position,res.ref_aa_position); 
        assert_eq!(mut_info.mut_aa_position,res.mut_aa_position); 
        assert_eq!(mut_info.ref_aa,res.ref_aa);
        assert_eq!(mut_info.mut_aa,res.mut_aa);
    }

    #[test]
    fn test_parse_amino_acid_seq_position()
    {
        let input_seq="32Q".to_string();
        match parse_amino_acid_seq_position(&input_seq)
        {
            Ok((pos,seq))=>
            {
                assert_eq!(pos,32u16);
                assert_eq!(seq,"Q".to_string());
            }
            Err(err_msg)=>
            {
                format!("Test failed, the code resturned an error message: {}",err_msg);
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_indepth()
    {
        let amino_acids_chars="ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>();
        for i in 0..100 // simulate different amino acids senario
        {
            let pos:u16=i;
            for j in 1..24 // simulate postential length
            {
                let random_seq=amino_acids_chars.iter().take(j).collect::<String>(); // get a test sequence
                let mut pos_dev=pos.to_string();
                pos_dev.push_str(&random_seq);
                match parse_amino_acid_seq_position(&pos_dev)
                {
                    Ok((res_pos,res_seq))=>
                    {
                        assert_eq!(pos,res_pos);
                        assert_eq!(random_seq,res_seq);
                    },
                    Err(err_msg)=>
                    {
                        println!("Test Case faile with the following error message: {}, input string was: {}",err_msg,pos_dev);
                    }
                }
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_bad_input_1()->Result<(),String>
    {
        let test_case="Test"; // here test case should fail because there is no position  
        match parse_amino_acid_seq_position(&test_case)
        {
            Ok((pos,seq))=>
            {
                println!("Test Case failed it returned the following results: {},{}",pos,seq);
                Err(format!("Test Case failed it returned the following results: {},{}",pos,seq))
            },
            Err(err_msg)=>
            {
                println!("Test Case faile with the following error message: input string was: {}",err_msg);
                Ok(())
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_bad_input_2()->Result<(),String>
    {
        let test_case=""; // here test case should fail because there is no position  
        match parse_amino_acid_seq_position(&test_case)
        {
            Ok((pos,seq))=>
            {
                println!("Test Case failed it returned the following results: {},{}",pos,seq);
                Err(format!("Test Case failed it returned the following results: {},{}",pos,seq))
            },
            Err(err_msg)=>
            {
                println!("Test Case faile with the following error message: input string was: {}",err_msg);
                Ok(())
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_bugged_input()->Result<(),String>
    {
        let test_case="-223QK"; // here test case should fail because there is no position  
        match parse_amino_acid_seq_position(&test_case)
        {
            Ok((pos,seq))=>
            {
                println!("Test Case failed it returned the following results: {},{}",pos,seq);
                Err(format!("Test Case failed it returned the following results: {},{}",pos,seq))
            },
            Err(err_msg)=>
            {
                println!("Test Case faile with the following error message: input string was: {}",err_msg);
                Ok(())
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_rare_input_1()
    {
        let input_seq="32*".to_string();
        match parse_amino_acid_seq_position(&input_seq)
        {
            Ok((pos,seq))=>
            {
                assert_eq!(pos,32u16);
                assert_eq!(seq,"*".to_string());
            }
            Err(err_msg)=>
            {
                println!("Test failed, the code resturned an error message: {}",err_msg);
            }
        }
    }
    #[test]
    fn test_parse_amino_acid_seq_position_rare_input_2()
    {
        let input_seq="32KMNOPQQQ*".to_string();
        match parse_amino_acid_seq_position(&input_seq)
        {
            Ok((pos,seq))=>
            {
                assert_eq!(pos,32u16);
                assert_eq!(seq,"KMNOPQQQ*".to_string());
            }
            Err(err_msg)=>
            {
                println!("Test failed, the code resturned an error message: {}",err_msg);
            }
        }
    }
    #[test]
    fn test_get_bit_mask(){}

    #[test]
    fn test_remove_leading_zeros_1()
    {
        let test_case1="3,4,0"; 
        let results=remove_leading_zeros(test_case1.to_string()); 
        assert_eq!(results,"3,4"); 
    }
    #[test]
    fn test_remove_leading_zeros_2()
    {
        let test_case1="3,4,0,1,0"; 
        let results=remove_leading_zeros(test_case1.to_string()); 
        assert_eq!(results,"3,4,0,1"); 
    }
    #[test]
    fn test_remove_leading_zeros_3()
    {
        let test_case1="0,0"; 
        let results=remove_leading_zeros(test_case1.to_string()); 
        assert_eq!(results,""); 
    }
    #[test]
    fn test_parse_fields1()
    {
        let test_case="0"; 
        let results=parse_fields(test_case.to_string());
        assert_eq!(results,"0$")
    }
    #[test]
    fn test_parse_fields2()
    {
        let test_case="6"; 
        let results=parse_fields(test_case.to_string());
        assert_eq!(results,"6$")
    }
    #[test]
    fn test_parse_fields3()
    {
        let test_case="6,3"; 
        let results=parse_fields(test_case.to_string());
        assert_eq!(results,"")
    }
    #[test]
    fn test_get_bit_mask1()
    {
        let test_case="0|0"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"");
    }
    #[test]
    fn test_get_bit_mask2()
    {
        let test_case="0|0:.:79,0:79:99:.:.:.:0"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"0$");
    }
    #[test]
    fn test_get_bit_mask3()
    {
        let test_case="0|0:.:37,0:37:72:.:.:.:0"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"0$");
    }
    #[test]
    fn test_get_bit_mask4()
    {
        let test_case="0|0:0"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"0$");
    }
    #[test]
    fn test_get_bit_mask5()
    {
        let test_case="0|1:0.541667:26,22:48:PASS:99:577,0,683:..:0.3336:2"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"2$");
    }
    #[test]
    fn test_get_bit_mask6()
    {
        let test_case="0|1:10"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"10$");
    } 
    #[test]
    fn test_get_bit_mask7()
    {
        let test_case="0|1:0.432432:16,21:37:PASS:99:634,0,417:..:0.1989:10922"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"10922$"); 
    }
    #[test]
    fn test_get_bit_mask8()
    {
        let test_case="1|1:.:4,87:91:99:3000,249,0:..:0.4777:15"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"15$"); 
    }
    #[test]
    fn test_get_bit_mask9()
    {
        let test_case="1|1:.:4,87:91:99:3000,249,0:..:0.4777:15,32,14"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"15,32,14"); 
    }
    #[test]
    fn test_get_bit_mask10()
    {
        let test_case="1|1:.:4,87:91:99:3000,249,0:..:0.4777:15,32,14,0,0,0"; 
        let results=get_bit_mask(&test_case.to_string());
        assert_eq!(results,"15,32,14"); 
    }
    #[test]
    fn test_get_types()
    {
        let test_case="*missense|ITPRID1|ENST00000409210|protein_coding|+|717C>717Y|31643796G>A".to_string(); 
        assert_eq!(*"*missense", *get_type(&test_case));
    }
}
