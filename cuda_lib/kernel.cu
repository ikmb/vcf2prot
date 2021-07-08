__global__
void apply_instruction_set_to_reference_seq_gpu(char* res_array, char* ref_stream, char* alt_Stream, 
    size_t* exc_code, size_t* start_pos, size_t* length, size_t* start_pos_res, unsinged int int num_task)
{
    int thread_index=(blockIdx.x * blockDim.x) + threadIdx.x; // get the thread index
	int step_size = blockDim.x + gridDim.x; // the step size

    for(int work_idx=thread_index; work_idx < num_task; work_idx +=step_size)
    { 
        if(exc_code[work_idx]==0) // A copy based instruction
        {
        // apply the copy based instruction for the whole length of the instruction
            for(int ins_count=0; ins_count!=length[work_idx];ins_count++)
            {
                res_array[start_pos_res[work_idx]+ins_count]=ref_stream[start_pos[work_idx]+ins_count]; 
            }
        }
        else // an insertion based instruction
        {
            for(int ins_count=0; ins_count!=instruction_length[work_idx];ins_count++)
            {
                res_array[start_pos_res[work_idx]+ins_count]=alt_Stream[start_pos[work_idx]+ins_count];  
            }
        }
    }

/** C wrapper */
extern "C"
{  
    /*@brief: A C wrapper for the CUDA kernel that perform Task execution 
    * @param ref_seq: A pointer to a char array the will hold the results 
    * @param ref_seq: A pointer to a char array holding the reference stream 
    * @param alt_stream: A pointer to a char array holding the alteration stream 
    * @param exc_code: A pointer to a size_t array holding the execution code results 
    * @param start_pos: A pointer to a size_t array holding the start position in the input stream, whether it is the ref or alternative
    * @param length: A pointer to a size_t array holding the length of the instruction 
    * @param start_pos_ref: A pointer to a size_t array holding the position in the resulting array  
    * @param num_task: the number of tasks which equal the length of the following arrays: ref_seq, ref_stream, alt_stream and exc_code
    * @param len_res_array: the length of the results array
    * @param len_ref_stream: the length of the refernce stream 
    * @param len_alt_stream: the length of the alt stream
    * @Notes: Error code meaning: 
    *       1. 0 => Sucess 
    *       2. 1 => GPU allocation failure 
    *       3. 2 => Failure with copying the data to the GPU 
    *       4. 3 => Launching the kernel failed 
    *       5. 4 => Kernel execution failed
    *       6. 5 => Copying the results array to the host failed 
    */
    int kernel_wrapper(char* res_array, char* ref_stream, char* alt_stream, size_t* exc_code,
        size_t* start_pos, size_t* length, size_t* start_pos_res, size_t num_task, 
        size_t len_res_array, size_t len_ref_stream,
        size_t len_alt_stream) 
    {
        // allocate arrays on the GPU
        //---------------------------
        // 1. creating pointers:
        //----------------------
        char* res_array_ptr; 
        char* ref_stream_ptr; 
        char* alt_stream_ptr; 
        size_t* exc_code_ptr; 
        size_t* start_pos_ptr; 
        size_t* length_ptr; 
        size_t* start_pos_res_ptr; 
        // 2. Perform allocations 
        //-----------------------
        if(cudaMalloc(&res_array_ptr,len_res_array*sizeof(char))!=cudaSuccess)return 1;
        if(cudaMalloc(&ref_stream_ptr,len_ref_stream*sizeof(char))!=cudaSuccess)return 1; 
        if(cudaMalloc(&alt_stream_ptr,len_alt_stream*sizeof(char))!=cudaSuccess)return 1; 
        if(cudaMalloc(&exc_code_ptr,num_task*sizeof(size_t))!=cudaSuccess)return 1; 
        if(cudaMalloc(&start_pos_ptr,num_task*sizeof(size_t))!=cudaSuccess)return 1; 
        if(cudaMalloc(&length_ptr,num_task*sizeof(size_t))!=cudaSuccess)return 1; 
        if(cudaMalloc(&start_pos_res_ptr,num_task*sizeof(size_t))!=cudaSuccess)return 1; 
        // 3. Copy the data to the GPU
        //---------------------------- 
        if(cudaMemcpy(ref_stream_ptr, ref_stream, len_ref_stream*sizeof(char),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        if(cudaMemcpy(alt_stream_ptr, alt_stream, len_alt_stream*sizeof(char),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        if(cudaMemcpy(exc_code_ptr, exc_code, num_task*sizeof(size_t),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        if(cudaMemcpy(start_pos_ptr, start_pos, num_task*sizeof(size_t),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        if(cudaMemcpy(length_ptr, length, num_task*sizeof(size_t),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        if(cudaMemcpy(start_pos_res_ptr, start_pos_r, num_task*sizeof(size_t),cudaMemcpyHostToDevice)!=cudaSuccess)return 2;
        // Launching the kernel 
        //----------------------
        unsigned int num_threads_per_block=1024; 
        unsigned int number_blocks=(num_task/num_threads_per_block) +1 ; 
        apply_instruction_set_to_reference_seq_gpu<<<number_blocks,num_threads_per_block>>>
        (
            res_array_ptr, ref_stream_ptr, alt_stream_ptr, exc_code_ptr, 
            start_pos_ptr, length_ptr, start_pos_res_ptr, num_task 
        );         
        if(cudaGetLastError()!=cudaSuccess)return 3;
        // Synchronize the calling code: 
        //------------------------------
        cudaErr_t err=cudaDeviceSynchronize(); 
        if (err!=cudaSuccess) return 4; 
        // Copying the from the GPU 
        //-------------------------
        if(cudaMemcpy(res_array, res_array_ptr, len_res_array*sizeof(char), cudaMemcpyDeviceToHost)!=cudaSuccess) return 5;
        // Relase allocated on the GPU 
        //----------------------------
        cudaFree(res_array_ptr); 
        cudaFree(ref_stream_ptr); 
        cudaFree(alt_stream_ptr); 
        cudaFree(exc_code_ptr); 
        cudaFree(start_pos_ptr); 
        cudaFree(length_ptr); 
        cudaFree(start_pos_res_ptr);
        // return a success state
        //------------------------
        return 0;         
    }
} 