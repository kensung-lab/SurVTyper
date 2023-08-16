# SurVTyper
A genotyper for structural variations using NGS paired-end sequencing data.

## Installation

Download and uncompress the latest release. From inside the folder, run
```
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## Running

SurVTyper needs a VCF file with the SVs to genotype, BAM/CRAM file, a (possibly empty) working directory and reference genome in FASTA format. 

The command for running it is:
```
python3 survtyper.py --threads N_THREADS VCF_FILE BAM_FILE WORKDIR REFERENCE
```

The genotyped VCF will be placed in WORKDIR/genotyped.vcf.gz.
