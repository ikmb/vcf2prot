use libc::{c_int, size_t};
/// A wrapper for C-linker and the executioner of the CUDA kernel  
extern "C" 
{
    pub fn kernel_wrapper(res_array:*mut u8, ref_stream: *const u8, alt_stream: *const u8, 
    exe_code: *const size_t, start_pos: *const size_t, length: *const size_t, start_pos_res: *const size_t,
    num_taks: size_t, len_res_array: size_t, len_ref_stream: size_t, len_alt_stream: size_t)->c_int;
} 
