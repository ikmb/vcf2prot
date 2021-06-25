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
./PPGG_rust -f input_mutation_info.vcf -r reference_sequences.fasta -o personalized_proteomes_dir -seiv  
```

### Compilation from source ### 

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

### Singularity Container ###

To be added later 

### Citation ###

TBD

### Contact ### 
For further questions, please feel free to open an issue here or send an email to the developers at h.elabd@ikmb.uni-kiel.de

### Funding ###
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’) 

![IKMB_LOGO](/media/RTG1743.png)
