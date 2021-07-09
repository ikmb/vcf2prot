![IKMB_LOGO](/media/IKMB_LOGO.png)

# PPGG: Personalized Proteome Generation using Graphical Processing Cards (GPUs) #

### Project Aim ### 

Accelerate the generation of personalized proteomes from a Variant calling format (VCF) file and a reference proteome using graphical processing units (GPUs).   

### Motivation ### 

To Be added later 

### Usage ### 

<p> Two mandatory inputs are needed by the tool, the first is the VCF containing the consequences calling and the second is a FASTA file containing reference sequences. </p>

#### NOTE ####  
<p> The program assumes the FASTA file to have the following structure </p>

```
>TRANS_ID
TRANS_SEQ_LINE1
TRANS_SEQ_LINE2 
>TRANS_ID
TRANS_SEQ_LINE1
.
.
.
```

<p> That is, the parser expects every char between > and '\n' to be the transcript name. Also, please make sure that ids used in the file are the same as in the VCF files. Otherwise, the program will not be able to functional properly. </p>

<p> Once the inputs have been defined the program can be used as follow: </p>

```
./PPGG_rust -f input_mutation_info.vcf -r reference_sequences.fasta -o personalized_proteomes_dir -sv  
```

## Compilation from source ## 

### CPU Version ###

#### **Note** ####
<p> Compiling the following code will be produce a CPU only version, that means that providing the code with panic if the GPU is specified as an engine, i.e. the parameter -g is set to gpu. </p>

<p> 1- Install Rust from https://www.rust-lang.org </p>

<p> 2- Clone the current repository </p>

```
git clone https://github.com/ikmb/ppg
```

<p> 3- Change the direction to ppg</p>

```
cd ppg
```

<p> 4- build the project </p>

```
cargo build 
```
<p> 5- Access the binary executable from the target directory </p> 

### GPU Version ###

#### **Note** ####

<p> The following GPU code is only compatible with CUDA and NVIDIA GPUs</p>

<p> 1- Install Rust from https://www.rust-lang.org </p>

<p> 2- Clone the current repository or Download the source code using the project Github page</p>

```
git clone https://github.com/ikmb/ppg
```

<p> 3- Change the direction to ppg</p>

```
cd ppg
```
<p> 4. Please make sure the following environmental variable are set CUDA_HOME and LD_LIBRARY_PATH, please set the value of these according to your system. </p>

<p> 5. Use any text editor and update the following information in the build script, *build.rs* which is located the at the root directory, the following the 8th</p>

´´´rust 
    println!("cargo:rustc-link-search=native=/opt/cuda/11.0/lib64/"); // 8th line in the current version
    println!("cargo:rustc-link-search=native=/path two cuda lib64 directory"); // 8th line in the updated version
´´´
<p> 4- build the project </p>

```
cargo build 
```
<p> 5- Access the binary executable from the target directory </p>

### Troubleshooting ###

#### Problem #### 

<p> error while loading shared libraries: libcudart.so.11.0: cannot open shared object file: No such file or directory </p>

#### solution #### 

<p> This problem will be encountered in case any of the two environmental variable, CUDA_HOME and LD_LIBRARY_PATH, are not defined or set. For a permanent solution please update your .bashrc to have these two variables exported.</p>

### Citation ###

TBD

### Contact ###
For further questions, please feel free to open an issue here or send an email to the developers at h.elabd@ikmb.uni-kiel.de

### Funding ###
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’)

![IKMB_LOGO](/media/RTG1743.png)