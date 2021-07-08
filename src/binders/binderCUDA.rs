use libc::{c_int, size_t};

#[link(name="kernel_wrapper", kind="static")]
extern "C" 
{
    pub fn kernel_wrapper(res_array:*mut char, ref_stream: *const char, alt_stream: *const char, 
    exe_code: *const size_t, start_pos: *const size_t, length: *const size_t, start_pos_res: *const size_t,
    num_taks: size_t, len_res_array: size_t, len_ref_stream: size_t, len_alt_stream: size_t)->c_int;
} 
