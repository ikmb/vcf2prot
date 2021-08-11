![IKMB_LOGO](/media/IKMB_LOGO.png)

# PPGG: Personalized Proteome Generation using Graphical Processing Cards (GPUs) #

## Project Aim ##

Accelerate the generation of personalized proteomes from a Variant calling format (VCF) file and a reference proteome using graphical processing units (GPUs).

## Execution Logic and Requirements ##

### Input Requirements ###

1. A reference fasta file containing transcript ids as sequence identifiers and the protein sequences of each transcript, for example, 

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

2. A VCF file containing the variants observed in the study population. The VCF file should be generated by  <a href= "https://academic.oup.com/bioinformatics/article/33/13/2037/3000373"> BCF/csq </a> as PPGG has been optimized to decode it's bit-mask and to parse it's consequence field. The file should also be phased and in a flat-VCF not BCF format.

#### Notes ####

1. The only exception is this is when the python wrapper is used which work directly with BCF tabix indexed files. 

2. You can decode a BCF file into a VCF using the following command:

```
bcftools view PATH_TO_BCF -O v -o PATH_TO_VCF
```

### Hardware Requirements ###

#### GPU version ####

<p> The GPU version of PPGG expects Nvidia-GPU to be accessible on the system during the development we utilized Tesla V100 SXM2 32GB. </p>

#### CPU version ###

<p> Expects a modern multi-core CPU with a big enough RAM to hold the whole file in memory during development a compute node with 512 GB of RAM and a twin intel Xeon CPU were used. </p>

### Software Requirements ###

The GPU version of the code can be compiled on a Linux-system with an available NVCC compliler and an Nvidia GPU.

The CPU version of the code can be compiled on a Linux and Mac OS system with Cargo.  

### Execution Logic ###

PPGG execution logic can be separate into the following main steps:

1. Reading and parsing the file where the file is read as a UTF-8 encoded string, patient names are extracted and records are filtered where only record with a supported protein coding effect are included into the next step. List of alterations supported by the current mutation is available in the file list_supported_alterations.tsv.

2. Once the VCF Records have been filtered, bit-masks are decoded and combined with the consequence mutation to generate a hashmap linking each patient to a collection of mutation observed in both of the patients haplotypes.

3. For each patient, mutations are grouped by the transcript id, i.e. all mutation occurring on a specific transcript are combined together.

4. For each collection of mutations, mutations are translated into instructions, at that stage mutations are checked for logical errors, e.g.  Mutational Engulfment, Where one mutation is a subset of another mutation, or Multiple annotations, where for the same position is annotated with more than one mutation. Also semantic-equivalence where two mutations are different at the genetic level but are equivalent at the protein level is taken place leading to a much smaller and a more consistence definition of alterations at the protein-level. In case any logical error was encountered, a waring message is printed to the standard output descriptor and the transcript is filtered out. Finally, instructions are interpreted and a simple representation for the sequence transcript is generated, internally, this is represented a vector of Tasks.  

5. After encoding each transcript into tasks, all transcripts are concatenated end-to-end to generate a vector of tasks describing the generation of all sequences in the haplotype.

6. Next, a backend engine is used to execute the tasks and generate the sequences for example, this engine can be a collection of CPU-threads or an execution stream on the GPU.

7. Finally, the results of the file are written to the Desk using a pool of writer-threads

### Usage ###

<p> Two mandatory inputs are needed by the tool, the first is the VCF containing the consequences calling and the second is a FASTA file containing reference sequences. </p>

#### Example ####

#### Clone the project ####

```bash
git clone https://github.com/ikmb/ppg
```

<p> Please note that git usually comes  pre-installed on most Mac OS and Linux systems. If git is not available at your system, you can install it from <a href= "https://git-scm.com/book/en/v2/Getting-Started-Installing-Git"> here </a>   </p>

#### Change directory to the ppg ####

```bash
cd ppg 
```

<p> Please notice that after calling git, a directory named ppg in the directory from which git has been called.</p>

##### Is ppg installed ? #####

To follow along, make sure the executable ppg has been installed in your system and is available on your PATH. Incase it is not installed, check the installation guideline below.

##### Export Env variables #####

<p> Let's  Inspect the GPU arrays, instruction's generation and the Task's arrays </p>

```bash
export DEBUG_GPU=TRUE
export INSPECT_TXP=TRUE
export INSPECT_INS_GEN=TRUE
```

<p> for more details about the meaning of the exported, check the Environment Variables section below </p>

##### Create a new directory to store the results #####

```bash
mkdir results 
```

##### Call PPG with the generated data #####

```bash
ppg -f examples/example_file.vcf -r examples/References_sequences.fasta -vs -g st -o results
```

<p> Where o flag determine the path to write the fasta file, the s guide the program to write stats and v for printing log statement. </p>

#### Environment Variables ####  

PPGG also utilize environmental variable heavily to customize its behavior, the list of environmental variable utilized by the  PPGG is shown below: 

1. DEBUG_GPU => Inspect the input arrays to the GPU are inspected for indexing error, incase of an indexing error the full input table is printed and idex
of the row with the first indexing error is also printed to the standard output descriptor.

2. DEBUG_CPU_EXEC => Inspect the vector of tasks provided to the input CPU execution engine, incase of an indexing error the full input table is printed and idex
of the row with the first indexing error is also printed to the standard output descriptor.

3. DEBUG_TXP="Transcript_ID" => This flag exports a transcript id that will be used for debugging, while the transcript id for transcript is being create different infos will be logged to the output descriptor.

4. INSPECT_TXP => If set, after each transcript is translated into instruction an inspection function will be called to check the correctness of translation, if the translation failed then the code will panic and error will be printed to the output descriptor.

5. INSPECT_INS_GEN => Inspect the translation process from mutations to instructions, as of version 0.1.3 two logical errors are inspected, first,
 multiple annotations, where more than one mutation are observed at the same position in the protein backbone, or through mutational overlap and engulfment where two mutations overlap in length, for example, insertion at position 60 with 7 amino acids and then a missense mutation at position 64.

6. PANIC_INSPECT_ERR => If set the code will panic if inspecting the translation from mutation to instruction failed. This is an override of the default behavior were an error message is generated and printed to the output stream.

## Compilation from source ##

### CPU Version ###

#### **Note** ####

<p> Compiling the following code will be produce a CPU only version, that means that providing the code with panic if the GPU is specified as an engine, i.e. the parameter -g is set to gpu. </p>

1. Install Rust from the <a href= "https://www.rust-lang.org "> official website </a>  

2. Clone the current repository

```bash
git clone https://github.com/ikmb/ppg
```

3. Change the direction to ppg

```bash
cd ppg
```

4. Change to the cpu-only branch

```bash
git checkout cpu-only
```

5. build the project

```bash
cargo build --release 
```

6. Access the binary executable from the target directory

```bash
cd target/release
./ppg -h # This print the help statement 
```

7. add the binary to your PATH

### GPU Version ###

<p> The following GPU code is only compatible with CUDA and NVIDIA GPUs</p>

1. Install Rust from the <a href= "https://www.rust-lang.org "> official website </a>  

2. Clone the current repository or Download the source code using the project Github page

```bash
git clone https://github.com/ikmb/ppg
```

3. Change the direction to ppg

```bash
cd ppg
```

4. Please make sure the following environmental variable are set CUDA_HOME and LD_LIBRARY_PATH, please set the value of these according to your system. 

5. Use any text editor and update the following information in the build script, *build.rs* which is located the at the root directory, the following the 8th

```rust
    println!("cargo:rustc-link-search=native=/opt/cuda/11.0/lib64/"); // 8th line in the current version
    println!("cargo:rustc-link-search=native=/path two cuda lib64 directory"); // 8th line in the updated version
```

6. build the project 

```bash
cargo build --release 
```

5. Access the binary executable from the target directory

```bash
cd target/release
./ppg -h # This print the help statement 
```

## Troubleshooting ##

### Problem ###

<p> error while loading shared libraries: libcudart.so.11.0: cannot open shared object file: No such file or directory </p>

### solution ###

<p> This problem will be encountered in case any of the two environmental variable, CUDA_HOME and LD_LIBRARY_PATH, are not defined or set. For a permanent solution please update your .bashrc to have these two variables exported.</p>


## Contact ##
For further questions, please feel free to open an issue here or send an email to the developers at h.elabd@ikmb.uni-kiel.de or through twitter @HeshamElAbd16

## Funding ##
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’)

![IKMB_LOGO](/media/RTG1743.png)