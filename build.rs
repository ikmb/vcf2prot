fn main() 
{
    cc::Build::new()
        .cuda(true)
        .flag("-cudart=shared")
        .file("cuda_lib/kernel.cu")
        .compile("libcuda_engine");
    println!("cargo:rustc-link-search=native=/opt/cuda/11.0/lib64/");
    println!("cargo:rustc-link-lib=dylib=cudart");
    
}