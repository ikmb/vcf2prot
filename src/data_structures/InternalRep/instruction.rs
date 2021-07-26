use crate::data_structures::{Constants, mutation_ds::*, vcf_ds::AltTranscript}; 
use rayon::vec;
use serde::{Deserialize, Serialize};
// the automatic derivatization of traits 
#[derive(Debug,Clone,Serialize,Deserialize,PartialEq)]
pub struct Instruction
{
    code:char,
    s_state:bool,
    pos:usize,
    len:usize,
    data:Vec<char>
}
impl Instruction
{
    /// create a new instruction where code is the instruction code, pos, is the position code,
    /// len is the length of the instruction and finally, data is the data associated with the 
    /// instruction, for example, the sequence of an insertion and/or a deletion. 
    /// The meaning of each code is as follow: 
    /// Missense ->	    M
    /// *Missense ->    N
    /// frameshift ->	F
    /// *frameshift	->  R
    /// stop_gained ->	G
    /// Stop_lost ->	L
    /// inframe_insertion ->    I	
    /// *inframe_insertion ->	J
    /// inframe_deletion -> 	D
    /// *inframe_deletion ->	C
    /// *missense&inframe_altering  ->	    K
    /// *frameshift&stop_retained ->	    Q
    /// *stop_gained&inframe_altering ->	A
    /// frameshift&stop_retained ->     	B
    /// inframe_deletion&stop_retained ->	P
    /// inframe_insertion&stop_retained	->  Z
    /// stop_gained&inframe_altering -> 	T
    /// stop_lost&frameshift ->             W
    /// missense&inframe_altering ->        Y
    /// start_lost&splice_region ->         U
    pub fn new(code:char,s_state:bool,pos:usize,len:usize,data:Vec<char>)->Self
    {
        Instruction{code, s_state, pos, len, data}
    }
    /// this is going to be the main translator of the language, it takes as input the mutation type 
    /// an returns an instruction Representing the interpreted code 
    pub fn from_mutation(mutation:&Mutation, vec_mut:&Vec<Mutation>)->Self
    {
        match &mutation.mut_type
        {
            MutationType::MisSense=>Instruction::interpret_missense(mutation,vec_mut), 
            MutationType::SMisSense=>Instruction::interpret_s_missense(mutation,vec_mut),
            MutationType::FrameShift=>Instruction::interpret_frameshift(mutation,vec_mut),
            MutationType::SFrameShift=>Instruction::interpret_s_frameshift(mutation,vec_mut),
            MutationType::InframeInsertion=>Instruction::interpret_inframe_insertion(mutation,vec_mut),
            MutationType::SInframeInsertion=>Instruction::interpret_s_inframe_insertion(mutation,vec_mut),
            MutationType::InframeDeletion=>Instruction::interpret_inframe_deletion(mutation, vec_mut),
            MutationType::SInframeDeletion=>Instruction::interpret_s_inframe_deletion(mutation, vec_mut),
            MutationType::StartLost=>Instruction::interpret_start_lost(mutation,vec_mut),
            MutationType::StopLost=>Instruction::interpret_stop_lost(mutation,vec_mut),
            MutationType::StopGained=>Instruction::interpret_stop_gained(mutation,vec_mut),
            MutationType::SStopGained=>Instruction::interpret_s_stop_gained(mutation,vec_mut), 
            MutationType::SMisSenseAndInframeAltering=>Instruction::interpret_s_missense_and_inframe_altering(mutation,vec_mut), 
            MutationType::SFrameShiftAndStopRetained=>Instruction::interpret_s_frameshift_and_stop_retained(mutation,vec_mut),
            MutationType::SStopGainedAndInframeAltering=>Instruction::interpret_s_stop_gained_and_inframe_altering(mutation,vec_mut), 
            MutationType::FrameShiftAndStopRetained=>Instruction::interpret_frameshift_and_stop_retained(mutation, vec_mut), 
            MutationType::InframeDeletionAndStopRetained=>Instruction::interpret_inframe_deletion_and_stop_retained(mutation,vec_mut),
            MutationType::InframeInsertionAndStopRetained=>Instruction::interpret_inframe_insertion_and_stop_retained(mutation),
            MutationType::StopGainedAndInframeAltering=>Instruction::interpret_stop_gained_and_inframe_altering(mutation,vec_mut),
            MutationType::StopLostAndFrameShift=>Instruction::interpret_stop_lost_and_frameshift(mutation,vec_mut), 
            MutationType::MissenseAndInframeAltering=>Instruction::interpret_missense_and_inframe_altering(mutation,vec_mut),
            MutationType::StartLostAndSpliceRegion=>Instruction::interpret_start_lost_and_splice_region(mutation,vec_mut),        
        }
    }
    /// return the code of the instruction
    pub fn get_code(&self)->char
    {
        self.code
    }
    /// return the position of the instruction in the reference code 
    pub fn get_position(&self)->usize
    {
        self.pos
    }
    /// return the length of the instruction  
    pub fn get_length(&self)->usize
    {
        self.len
    }
    /// return a copy of the associated data 
    pub fn get_data(&self)->Vec<char>
    {
        self.data.clone()
    }
    pub fn get_s_state(&self)->bool
    {
        self.s_state
    }
    pub fn update_code(&mut self, code:char)
    {
        self.code=code;
    }
    pub fn update_s_state(&mut self,s_state:bool)
    {
        self.s_state=s_state; 
    }
    pub fn update_start_pos(&mut self, start_pos:usize)
    {
        self.pos=start_pos; 
    }
    /// add parser function that generate an instruction from a specific mutation type 
    fn generate_phi_instruction()->Self
    {
        let code='E'; 
        let len=0;
        let pos=0;
        let data=Vec::new(); 
        let s_state=false;
        Instruction{code, s_state, pos, len, data}
    }
    fn interpret_missense(mutation:&Mutation,_vec_alt_tram:&Vec<Mutation>)->Self
    {   
        let code='M'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid 
        let data= match &mutation.mut_info.mut_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data
            }
            MutatedString::NotSeq =>panic!("Something went wrong, interpreting: {:#?}, failed",&mutation)
        }; 
        let len=1;
        let s_state=false;
        Instruction{code, s_state, pos,len, data}
    }
    fn interpret_s_missense(mutation:&Mutation,vec_alt_tram:&Vec<Mutation>)->Self
    {
        match Instruction::validate_s_state(mutation,vec_alt_tram)
        {
            true=>
            {
                let mut n_inst=Instruction::interpret_missense(mutation,vec_alt_tram);
                let pos=
                n_inst.update_code('N'); 
                n_inst.update_s_state(true); 
                n_inst
            },
            false=>
            {
                Instruction::generate_phi_instruction()
            }
        }
    }
    fn interpret_inframe_insertion(mutation:&Mutation,_vec_alt_tram:&Vec<Mutation>)->Self
    {
        let code='I'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid
        let data= match &mutation.mut_info.mut_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data
            }
            MutatedString::NotSeq =>panic!("Something went wrong, interpreting: {:#?}, failed",&mutation)
        }; 
        let len=data.len();
        let s_state=false;
        Instruction{code, s_state, pos,len, data}
    }
    fn interpret_s_inframe_insertion(mutation:&Mutation,vec_alt_tram:&Vec<Mutation>)->Self
    {
        match Instruction::validate_s_state(mutation,vec_alt_tram)
        {
            true=>
            {
                let mut n_inst=Instruction::interpret_inframe_insertion(mutation,vec_alt_tram);
                n_inst.update_code('J'); 
                n_inst.update_s_state(true); 
                n_inst
            }
            false=>
            {
                Instruction::generate_phi_instruction()
            }
        }
    }
    fn interpret_inframe_deletion(mutation:&Mutation,_vec_mut:&Vec<Mutation>)->Self
    {
        let code='D'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid
        let len = match &mutation.mut_info.ref_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>().len(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data.len()
            }
            MutatedString::NotSeq => panic!("Something went wrong, interpreting: {:#?}, failed",&mutation)
        }; 
        let data=match &mutation.mut_info.mut_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data
            }
            MutatedString::NotSeq => panic!("Something went wrong, interpreting: {:#?}, failed",&mutation)
        }; 
        // the length of deletion is 1.
        let s_state=false;
        Instruction::new(code, s_state, pos, len - data.len(), data)
    }
    fn interpret_s_inframe_deletion(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        match Instruction::validate_s_state(mutation,vec_mut)
        {
            true=>
            {
                let mut n_inst=Instruction::interpret_inframe_deletion(mutation,vec_mut);
                n_inst.update_code('C'); 
                n_inst.update_s_state(true); 
                n_inst
            },
            false=>
            {
                Instruction::generate_phi_instruction()
            }
        }
    }
    fn interpret_frameshift(mutation:&Mutation, _vec_mut:&Vec<Mutation>)->Self
    {
        let code='F'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid
        let data= match &mutation.mut_info.mut_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data
            }
            MutatedString::NotSeq => return Instruction::generate_phi_instruction()
        }; 
        let len=data.len();// Becuase we have the first amino acid in the mutated sequences already, for example, 115SL>115S
        // the length of deletion is 1.
        let s_state=false;
        Instruction{code, s_state, pos, len, data}
    }
    fn interpret_s_frameshift(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        match Instruction::validate_s_state(mutation,vec_mut)
        {
            true=>
            {
                let mut n_inst=Instruction::interpret_frameshift(mutation,vec_mut);
                n_inst.update_code('R'); 
                n_inst.update_s_state(true); 
                n_inst
            },
            false =>
            {
                Instruction::generate_phi_instruction()
            }
        }
    }
    fn interpret_stop_gained(mutation:&Mutation, vec_mut:&Vec<Mutation>)->Self
    {
        let code='G'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid 
        let len=0;
        let data=Vec::new(); 
        let s_state=false;
        Instruction{code, s_state, pos, len, data}
    }
    fn interpret_s_stop_gained(mutation:&Mutation, vec_mut:&Vec<Mutation>)->Self
    {
        match Instruction::validate_s_state(mutation,vec_mut)
        {
            true=>
            {
                let mut n_inst=Instruction::interpret_stop_gained(mutation,vec_mut);
                n_inst.update_code('X'); 
                n_inst.update_s_state(true); 
                n_inst
            },
            false=>
            {
                Instruction::generate_phi_instruction()
            }
        }
    }
    fn interpret_stop_lost(mutation:&Mutation, _vec_mut:&Vec<Mutation>)->Self
    {
        let code='L'; 
        let pos=mutation.mut_info.mut_aa_position as usize; // the position of the altered amino acid 
        let data= match &mutation.mut_info.mut_aa
        {
            MutatedString::Sequence(seq_str)=>seq_str.chars().collect::<Vec<char>>(),
            MutatedString::EndSequence(seq_str)=>
            {
                let mut data=seq_str.chars().collect::<Vec<char>>();
                data.remove(data.len()-1);   
                data
            }
            MutatedString::NotSeq => panic!("Something went wrong, interpreting: {:#?}, failed",&mutation)
        }; 
        let len=data.len();
        let s_state=false;
        Instruction{code, s_state, pos, len, data}
    }
    fn interpret_start_lost(_mutation:&Mutation, _vec_mut:&Vec<Mutation>)->Self
    {
        let code='0'; 
        let len=0;
        let pos=0;
        let data=Vec::new(); 
        let s_state=false;
        Instruction{code, s_state, pos, len, data}
    }
    fn interpret_s_missense_and_inframe_altering(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_s_frameshift(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('K'); n_inst}
        }
    }
    fn interpret_s_frameshift_and_stop_retained(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_s_frameshift(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('Q'); n_inst}
        }
    }
    fn interpret_s_stop_gained_and_inframe_altering(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_s_stop_gained(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('A'); n_inst}
        }
    }
    fn interpret_frameshift_and_stop_retained(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_frameshift(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('B'); n_inst}
        }
    }
    fn interpret_inframe_deletion_and_stop_retained(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_stop_gained(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('P'); n_inst}
        }
    }
    fn interpret_inframe_insertion_and_stop_retained(mutation:&Mutation)->Self
    {
        let mut n_inst=Instruction::generate_phi_instruction();
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>
            {
                n_inst.update_code('Z');
                n_inst.update_start_pos(mutation.mut_info.mut_aa_position as usize);
                n_inst
            }
        }
    }
    fn interpret_stop_gained_and_inframe_altering(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_stop_gained(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('T'); n_inst}
        }
    }
    fn interpret_stop_lost_and_frameshift(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_stop_lost(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('W'); n_inst}
        }
    }
    fn interpret_missense_and_inframe_altering(mutation:&Mutation,vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_frameshift(mutation,vec_mut);
        match n_inst.get_code()
        {
            'E'=>n_inst,
            _=>{n_inst.update_code('Y'); n_inst}
        }
    }
    fn interpret_start_lost_and_splice_region(mutation:&Mutation, vec_mut:&Vec<Mutation>)->Self
    {
        let mut n_inst=Instruction::interpret_start_lost(mutation,vec_mut);
        n_inst.update_code('U');
        n_inst
    }
    fn validate_s_state(mutation:&Mutation,vec_alt_tram:&Vec<Mutation>)->bool
    {
        let index=vec_alt_tram.iter()
                .position(|elem|elem==mutation).unwrap(); 
        let mut state=true; 
        for mutation in vec_alt_tram[..index].iter()
        { 
            if mutation.mut_type == MutationType::StopGained || mutation.mut_type == MutationType::FrameShift
            {state=false; break;}
            
        }
        state
    }
}
#[cfg(test)]
pub mod test_instructions
{
    use super::*; 
    #[test]
    pub fn test_interpretion_of_missense_mut()
    {
        let test_case=vec!["missense".to_string(),"ENST00000484547".to_string(), "32Q>32R".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=Instruction::interpret_missense(&test_mutation); 
        assert_eq!(ins.get_code(),'M'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),31); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
        assert_eq!(ins.get_data()[0],'R'); 
    }
    #[test]
    pub fn test_interpretion_of_missense_mut2()
    {
        let test_case=vec!["*missense".to_string(),"ENST00000484547".to_string(), "32Q>32R".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=Instruction::interpret_s_missense(&test_mutation); 
        assert_eq!(ins.get_code(),'N'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),31); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
        assert_eq!(ins.get_data()[0],'R'); 
    }
    #[test]
    #[should_panic]
    pub  fn test_interpretion_of_missense_mut3()
    {
        let test_case=vec!["*missense".to_string(),"ENST00000484547".to_string(), "32Q>32*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=Instruction::interpret_s_missense(&test_mutation); 
        assert_eq!(ins.get_code(),'N'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),31); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
        assert_eq!(ins.get_data()[0],'R'); 
    }
    #[test]
    #[should_panic]
    pub  fn test_interpretion_of_missense_mut4()
    {
        let test_case=vec!["*missense".to_string(),"ENST00000484547".to_string(), "3200>32M".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        let ins=Instruction::interpret_missense(&test_mutation); 
        assert_eq!(ins.get_code(),'M'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),31); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
        assert_eq!(ins.get_data()[0],'M'); 
    }
    #[test]
    pub  fn test_interpretion_of_inframe_insertion_mut1()
    {
        // 125Y>125YRR   
        let test_case=vec!["inframe_insertion".to_string(),"ENST00000484547".to_string(), "125Y>125YRR".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);
        let ins=Instruction::interpret_inframe_insertion(&test_mutation); 
        assert_eq!(ins.get_code(),'I'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),124); 
        assert_eq!(ins.get_length(),3); 
        assert_eq!(ins.get_data().len(),3); 
        assert_eq!(ins.get_data()[0],'Y'); 
        assert_eq!(ins.get_data()[1],'R'); 
        assert_eq!(ins.get_data()[2],'R'); 
    }
    #[test]
    pub  fn test_interpretion_of_s_inframe_insertion_mut2()
    {
        // 125Y>125YRR   
        let test_case=vec!["*inframe_insertion".to_string(),"ENST00000484547".to_string(), "125Y>125YRR".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);
        let ins=Instruction::interpret_s_inframe_insertion(&test_mutation); 
        assert_eq!(ins.get_code(),'J'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),124); 
        assert_eq!(ins.get_length(),3); 
        assert_eq!(ins.get_data().len(),3); 
        assert_eq!(ins.get_data()[0],'Y'); 
        assert_eq!(ins.get_data()[1],'R'); 
        assert_eq!(ins.get_data()[2],'R'); 
    }
    #[test]
    pub fn test_interpretion_of_inframe_deletion()
    {
        let test_case=vec!["inframe_deletion".to_string(),"ENST00000506382".to_string(), "115SL>115S".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_inframe_deletion(&test_mutation); 
        println!("{:#?}",&ins);  
        assert_eq!(ins.get_code(),'D'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),114); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
    }
    #[test]
    pub fn test_interpretion_of_s_inframe_deletion()
    {
        let test_case=vec!["*inframe_deletion".to_string(),"ENST00000506382".to_string(), "115SL>115S".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_inframe_deletion(&test_mutation);
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'C'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),114); 
        assert_eq!(ins.get_length(),1); 
        assert_eq!(ins.get_data().len(),1); 
    }
    #[test]
    pub fn test_interpretion_of_frameshift()
    {
        let test_case=vec!["frameshift".to_string(),"ENST00000510017".to_string(), "40VGLHFWTM*>40VDSTFGQC".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_frameshift(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'F'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),39); 
        assert_eq!(ins.get_length(),8); 
        assert_eq!(*ins.get_data(),['V','D','S','T','F','G','Q','C']); 
    }
    #[test]
    pub fn test_interpretion_of_s_frameshift()
    {
        let test_case=vec!["*frameshift".to_string(),"ENST00000510017".to_string(), "40VGLHFWTM*>40VDSTFGQC".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_frameshift(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'R'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),39); 
        assert_eq!(ins.get_length(),8); 
        assert_eq!(*ins.get_data(),['V','D','S','T','F','G','Q','C']); 
    }
    #[test]
    pub fn test_interpretion_of_stop_gained()
    {
        let test_case=vec!["stop_gained".to_string(),"ENST00000313766".to_string(), "217E>217*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_stop_gained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'G'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),216); 
        assert_eq!(ins.get_length(),0); 
        assert_eq!(ins.get_data().len(),0); 
    }
    #[test]
    pub fn test_interpretion_of_stop_lost()
    {
        let test_case=vec!["stop_lost".to_string(),"ENST00000650310".to_string(), "489*>489S".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_stop_lost(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'L'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),488); 
        assert_eq!(ins.get_length(),1);
        assert_eq!(ins.get_data().len(),1);
        assert_eq!(*ins.get_data(),['S']);
    }
    #[test]
    pub fn test_interpretion_of_start_lost()
    {
        let test_case=vec!["start_lost".to_string(),"ENST00000275358".to_string(), "1M>1K".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_start_lost(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'0'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),0); 
        assert_eq!(ins.get_length(),0);
        assert_eq!(ins.get_data().len(),0);
    }
    #[test]
    pub fn test_interpretion_of_s_stop_gained()
    {
        let test_case=vec!["*stop_gained".to_string(),"ENST00000313766".to_string(), "217E>217*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_stop_gained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'X'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),216); 
        assert_eq!(ins.get_length(),0); 
        assert_eq!(ins.get_data().len(),0); 
    }
    #[test]
    pub fn test_interpret_s_missense_and_inframe_altering()
    {
        let test_case=vec!["*missense&inframe_altering".to_string(),"ENST00000326303".to_string(), "188LAY>188LQS".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_missense_and_inframe_altering(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'K'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),187); 
        assert_eq!(ins.get_length(),3); 
        assert_eq!(ins.get_data().len(),3); 
        assert_eq!(*ins.get_data(),['L','Q','S']); 
    }
    #[test]
    pub fn test_interpret_s_frameshift_and_stop_retained()
    {
        let test_case=vec!["*frameshift&stop_retained".to_string(),"ENST00000438700".to_string(), "308GSLGMGQLLLRAKAMRLLYYLKTEDPEYDVQSKQWLTHLLDQFTNIKNILALKKIEVVHFTSLSRQLEFEATSVTVIPVFHLAYILIILFAVTSCFRFDCIRNKMCVAAFGVISAFLAVVSGFGLLLHIGVPFVIIVANSPFLILGVGVDDMFIMISAWHKTNLADDIRERMSNVYSKAAVSITITTITNILALYTGIMSSFRSVQCFCIYTGMTLLFCYFYNITCFGAFMALDGKREVVCLCWLKKADPKWPSFKKFCCFPFGSVPDEHGTDIHPISLFFRDYFGPFLTRSESKYFVVFIYVLYIISSIYGCFHVQEGLDLRNLASDDSYITPYFNVEENYFSDYGPRVMVIVTKKVDYWDKDVRQKLENCTKIFEKNVYVDKNLTEFWLDAYVQYLKGNSQDPNEKNTFMNNIPDFLSNFPNFQHDINISSSNEIISSRGFIQTTDVSSSAKKKILLF*>308GQPRNGPVTPAGQSHAAAVLPEDRGP*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_frameshift_and_stop_retained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'Q'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),307); 
        assert_eq!(*ins.get_data(),[
            'G','Q','P','R','N','G','P','V','T',
            'P','A','G','Q','S','H','A','A','A','V','L','P','E','D','R','G','P',]); 
    }
    #[test]
    pub fn test_interpret_s_stop_gained_and_inframe_altering()
    {
        let test_case=vec!["*stop_gained&inframe_altering".to_string(),"ENST00000275358".to_string(), "1273KEEDDKNAQGRKRHVKPTSGNAFTICKYPCGKSRECVAPNICKCKPGYIGSNCQTALCDPDCKNHGKCIKPNICQCLPGHGGATCDEEHCNPPCQHGGTCLAGNLCTCPYGFVGPRCETMVCNRHCENGGQCLTPDICQC>1273".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_s_stop_gained_and_inframe_altering(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'A'); 
        assert_eq!(ins.get_s_state(),true); 
        assert_eq!(ins.get_position(),1272); 

    }
    #[test]
    pub fn test_interpret_frameshift_and_stop_retained()
    {
        let test_case=vec!["frameshift&stop_retained".to_string(),"ENST00000381329".to_string(), "65IEREFENLYIENLELRREIDTLNERLAAEGQAIDGAELSKGQLKTKASHSTSQLSQKLKTTYKASTSKIVSSFKTTTSRAACQLVKEYIGHRDGIWDVSVAKTQPVVLGTASADHTALLWSIETGKCLVKYAGHVGSVNSIKFHPSEQLALTASGDQTAHIWRYAVQLPTPQPVADTSVSTFPYL*>65IENLKTFISKT*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_frameshift_and_stop_retained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'B'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),64); 
        assert_eq!(*ins.get_data(),['I','E','N','L','K','T','F','I','S','K','T']);
    }
    #[test]
    pub fn test_inframe_deletion_and_stop_retained()
    {
        let test_case=vec!["frameshift&stop_retained".to_string(),"ENST00000381329".to_string(), "733S*>733*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_inframe_deletion_and_stop_retained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'P'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_length(),0);
        assert_eq!(ins.get_position(),732); 
    }
    #[test]
    pub fn test_interpret_inframe_insertion_and_stop_retained()
    {
        let test_case=vec!["inframe_insertion&stop_retained".to_string(),"ENST00000551483".to_string(), "192*>192*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_inframe_insertion_and_stop_retained(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'Z'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_length(),0);
        assert_eq!(ins.get_position(),191); 
    }
    #[test]
    pub fn test_stop_gained_inframe_altering()
    {
        let test_case=vec!["stop_gained&inframe_altering".to_string(),"ENST00000328942".to_string(), "22LESVQCWIGIPFCAIYLIAMIGNSLLLSIIKSERSLHEPLYIFLGMLGATDIALASSIMPKMLGIFWFNVPEIYFDSCLLQMWFIHTLQGIESGILVAMALDRYVAICYPLRHANIFTHQLVIQIGTMVVLRAAILVAPCLVLIKCRFQFYHTTVISHSYCEHMAIVKLAAANVQVNKIYGLFVAFTVAGFDLTFITLSYIQIFITVFRLPQKEARFKAFNTCIAHICVFLQFYLLAFFSFFTHRFGS>22*".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_stop_gained_and_inframe_altering(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'T'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),21); 
    }
    #[test]
    pub fn test_stop_lost_and_frameshift()
    {
        let test_case=vec!["stop_lost&frameshift".to_string(),"ENST00000398786".to_string(), "134*>134N".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_stop_lost_and_frameshift(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'W'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),133); 
        assert_eq!(*ins.get_data(),['N']);  
    }
    #[test]
    pub fn test_interpret_start_lost_and_splice_region()
    {
        let test_case=vec!["start_lost&splice_region".to_string(),"ENST00000375110".to_string(), "1M>1I".to_string()];
        let test_mutation=Mutation::new(test_case).unwrap();
        println!("{:#?}",&test_mutation);  
        let ins=Instruction::interpret_start_lost_and_splice_region(&test_mutation); 
        println!("{:#?}",&ins); 
        assert_eq!(ins.get_code(),'U'); 
        assert_eq!(ins.get_s_state(),false); 
        assert_eq!(ins.get_position(),0); 
    }
}