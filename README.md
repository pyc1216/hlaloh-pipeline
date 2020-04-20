# For test

cd test
snakemake -s ../hla_loh.docker.smk loh_all -np



# NOTE
## data/
### abc_complete.fasta
This is from polysolver(conda or docker build)

### hla.dat
This is from LOHHLA(https://github.com/slagtermaarten/LOHHLA.git) /data/hla.dat

## images/
### polysolver.tar.gz
`docker pull sachet/polysolver:v4`
`docker save -o images/polysolver.tar.gz sachet/polysolver:v4`

### lohhla.tar.gz
`docker pull mleventhal/lohhla:latest`
`docker save -o images/lohhla.tar.gz mleventhal/lohhla:latest`
