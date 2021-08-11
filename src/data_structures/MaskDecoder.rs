use crate::data_structures::Constants; 

/// The Bitmask struct is an API for creating and handling bit-masks fields, used to index into the Consequences strings
/// and extract the consequence of the genomic alteration at a specific genomic location. 
///
#[derive(Debug,Clone)]
pub struct BitMask
{
    pub bitmask_elements:Option::<Vec<u32>>
}
impl BitMask
{
    /// Construct a new bit-mask instance from a string containing the bit-mask, the input string has been provided by the function
    /// get_bit_mask defined at the functions::text_parser module  
    /// ## Example
    ///``` 
    /// use ppg_rust::data_structures::MaskDecoder::BitMask;
    /// use ppg_rust::functions::text_parser; 
    /// let mut test_case="0|1:0.432432:16,21:37:PASS:99:634,0,417:..:0.1989:10922"; 
    /// let mut results=text_parser::get_bit_mask(&test_case.to_string());
    /// assert_eq!(results,"10922$"); 
    /// match BitMask::from_string(&mut results).bitmask_elements
    /// {
    ///    Some(vec)=>
    ///    {
    ///        assert_eq!(vec.len(),1);
    ///        assert_eq!(vec[0],10922);
    ///    },
    ///    None=>()
    /// }
    ///```
    pub fn from_string(input_string:&mut String)->Self
    {
        if *input_string==Constants::DEF_CONSEQ || input_string=="0$"
        {
            return BitMask{bitmask_elements:None}; 
        }
        if input_string.ends_with("$")
        {
            let input_string= input_string.strip_suffix("$").unwrap(); 
            let bitmask_vec:Vec<u32>=vec![input_string.parse::<u32>().unwrap()];
            return BitMask{bitmask_elements:Some(bitmask_vec)}
        }
        else
        {
            let bitmask_vec=input_string.split(",")
                                        .map(|elem| elem.parse::<u32>().unwrap())
                                        .collect::<Vec<u32>>();
            return BitMask{bitmask_elements:Some(bitmask_vec)}   
        }
    }
    /// Parse the u32 integer in the bitmask set and return a tuple of two vectors, the first vector contain the 
    /// index of CSQ observed in the first haplotype and the second contains the CSQ observed in the second haplotype.
    /// ## Example
    ///``` 
    /// use ppg_rust::data_structures::MaskDecoder::BitMask; 
    /// let mut test_case="3,3,3,3".to_string();
    /// let mut test_bitmask=BitMask::from_string(&mut test_case);
    /// match test_bitmask.get_indices()
    /// {
    ///    None=>(), 
    ///    Some((vec_h1,vec_h2))=>
    ///    {
    ///        assert_eq!(vec_h1.len(),4);
    ///        assert_eq!(vec_h2.len(),4);      
    ///        assert_eq!(vec_h1[0],0);
    ///        assert_eq!(vec_h1[1],15);
    ///        assert_eq!(vec_h1[2],30);
    ///        assert_eq!(vec_h1[3],45);
    ///        assert_eq!(vec_h2[0],0);
    ///        assert_eq!(vec_h2[1],15);
    ///        assert_eq!(vec_h1[2],30);
    ///        assert_eq!(vec_h1[3],45);
    ///    }
    /// }
    ///```
    pub fn get_indices(&mut self)->Option<(Vec<usize>,Vec<usize>)>
    {
        match &mut self.bitmask_elements
        {
            None=>None,
            Some(vec)=>
            {
                if vec.len()==1
                {
                    return Some(BitMask::parse_single_field(vec[0]));
                }
                else 
                {
                    return Some(BitMask::parse_concat_values(vec)); 
                }
            }
        }
    }
    fn parse_single_field(mut bitmask:u32)->(Vec<usize>,Vec<usize>)
    {
        let mut haplotype_one=Vec::with_capacity(16);
        let mut haplotype_two=Vec::with_capacity(16);
        let mut haplo1;
        let mut haplo2;  
        let mut index=0;
        // loop over all bits in the bitmask          
        while bitmask!=0
        {
            // decode the bit_mask 
            haplo1=bitmask;
            haplo2=bitmask>>1;
            if haplo1&1 == 1
            {
                haplotype_one.push(index);
            }
            if haplo2&1 == 1
            {
                haplotype_two.push(index);
            }
            // update the while loop element 
            bitmask=bitmask>>2;
            index+=1; 
        }
        (haplotype_one,haplotype_two)
    }
    fn parse_concat_values(bitmask_vec:&mut Vec<u32>)->(Vec<usize>,Vec<usize>)
    {
        let mut haplotype_one=Vec::with_capacity(16*bitmask_vec.len());
        let mut haplotype_two=Vec::with_capacity(16*bitmask_vec.len());
        let mut haplo1;
        let mut haplo2;  
        let mut index;
        let mut index_fields=0; 
        for bitmask in bitmask_vec
        {
            index=0;
            while *bitmask!=0
            {
                // decode the bit_mask 
                haplo1=*bitmask;
                haplo2=*bitmask>>1;
                if haplo1&1==1
                {
                    haplotype_one.push(index_fields+index);
                }
                if haplo2&1==1
                {
                    haplotype_two.push(index_fields+index);
                }
                // updat the while loop element 
                *bitmask=*bitmask>>2;
                index+=1;
            }
            index_fields+=15; 
        }
        (haplotype_one,haplotype_two)
    }
}

#[cfg(test)]
mod test_bitmask_class
{
    use super::*; 

    #[test]
    fn test_from_string1()->Result<(),()>
    {
        let mut test_case="".to_string(); 
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(_)=>Err(()),
            None=>Ok(())
        }
    }
    #[test]
    fn test_from_string2()->Result<(),()>
    {
        let mut test_case="2$".to_string();
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>
            {
                assert_eq!(vec.len(),1);
                assert_eq!(vec[0],2);
                Ok(())
            },
            None=>Err(())
        }

    }
    #[test]
    fn test_from_string3()->Result<(),()>
    {
        let mut test_case="0$".to_string();
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>Err(()),
            None=>Ok(())
        }
    }
    #[test]
    fn test_from_string4()->Result<(),()>
    {
        let mut test_case="1024$".to_string();
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>
            {
                assert_eq!(vec.len(),1);
                assert_eq!(vec[0],1024);
                Ok(())
            },
            None=>Err(())
        }
    }
    #[test]
    fn test_from_string5()->Result<(),()>
    {
        let mut test_case="1024,2048,4096".to_string();
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>
            {
                assert_eq!(vec.len(),3);
                assert_eq!(vec[0],1024);
                assert_eq!(vec[1],2048);
                assert_eq!(vec[2],4096);
                Ok(())
            },
            None=>Err(())
        }
    }
    #[test]
    fn test_from_string7()->Result<(),()>
    {
        let mut test_case="1024,0,4096".to_string();
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>
            {
                assert_eq!(vec.len(),3);
                assert_eq!(vec[0],1024);
                assert_eq!(vec[1],0);
                assert_eq!(vec[2],4096);
                Ok(())
            },
            None=>Err(())
        }
    }
    #[test]
    #[should_panic] 
    fn test_from_string8()->()
    {
        let mut test_case="-1,0,4096".to_string(); // if there is a negative number the function will panic !!! 
        match BitMask::from_string(&mut test_case).bitmask_elements
        {
            Some(vec)=>
            {
                assert_eq!(vec.len(),3);
                assert_eq!(vec[1],0);
                assert_eq!(vec[2],4096);
            },
            None=>()
        }
        ()
    }
    #[test]
    fn test_get_indicies1()->Result<(),()>
    {
        let mut test_case="0$".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Ok(()), 
            Some(vecs)=>Err(())
        }
    }
    #[test]
    fn test_get_indicies2()->Result<(),()>
    {
        let mut test_case="1$".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),1);
                assert_eq!(vec_h2.len(),0);
                assert_eq!(vec_h1[0],0);
                Ok(())
            }
        }
    }
    #[test]
    fn test_get_indicies3()->Result<(),()>
    {
        let mut test_case="3$".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),1);
                assert_eq!(vec_h2.len(),1);
                assert_eq!(vec_h1[0],0);
                assert_eq!(vec_h2[0],0);
                Ok(())
            }
        }
    }
    #[test]
    fn test_get_indicies4()->Result<(),()>
    {
        let mut test_case="1024$".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),1);
                assert_eq!(vec_h2.len(),0);      
                assert_eq!(vec_h1[0],5);
                Ok(())
            }
        }
    }
    #[test]
    fn test_get_indicies5()->Result<(),()>
    {
        let mut test_case="1,1".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),2);
                assert_eq!(vec_h2.len(),0);      
                assert_eq!(vec_h1[0],0);
                assert_eq!(vec_h1[1],15);
                Ok(())
            }
        }
    }
    #[test]
    fn test_get_indicies6()->Result<(),()>
    {
        let mut test_case="3,3".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),2);
                assert_eq!(vec_h2.len(),2);      
                assert_eq!(vec_h1[0],0);
                assert_eq!(vec_h1[1],15);
                assert_eq!(vec_h2[0],0);
                assert_eq!(vec_h2[1],15);
                Ok(())
            }
        }
    }
    #[test]
    fn test_get_indicies7()->Result<(),()>
    {
        let mut test_case="3,3,3,3".to_string();
        let mut test_bitmask=BitMask::from_string(&mut test_case);
        match test_bitmask.get_indices()
        {
            None=>Err(()), 
            Some((vec_h1,vec_h2))=>
            {
                assert_eq!(vec_h1.len(),4);
                assert_eq!(vec_h2.len(),4);      
                assert_eq!(vec_h1[0],0);
                assert_eq!(vec_h1[1],15);
                assert_eq!(vec_h1[2],30);
                assert_eq!(vec_h1[3],45);
                assert_eq!(vec_h2[0],0);
                assert_eq!(vec_h2[1],15);
                assert_eq!(vec_h1[2],30);
                assert_eq!(vec_h1[3],45);
                Ok(())
            }
        }
    }
}











