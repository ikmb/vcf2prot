# AUTHOR: Hesham ELAbd
# BREIF: Building the docker container for VCF2PROT 
# DATE: 15.07.2022

# Define the image root
FROM rust:1.61
# add the label and the version of the image 
LABEL Name=vcf2prot Version=0.1.4
# Copy the source code and the dependencies 
COPY ./src ./src
COPY Cargo.toml ./
# build the project with 
RUN cargo build --release 
# Run VCF2Prot
ENTRYPOINT ["./target/release/vcf2prot"]
