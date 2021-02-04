# fitDNM
### Evin Padhi
### 1/30/21
### Washington University in St. Louis
### School of Medicine, Department of Genetics
### Laboratory of Dr. Tychele Turner
### Snakefile for running fitDNM and all preprocessing steps
### Contact: evin.padhi@wustl.edu


### Purpose:
- To take in any genomic coordiantes, extract de novo variants from the callset and get scores every nucleotide and subsequently run through fitDNM that has been adapted to work in noncoding regions
- The goal is to be able to run fitDNM genome-wide

### How to run:
1. Download CADD scores and check md5sums
2. Build dockerfile and push
3. Modify .json  `fitDNM_genome_wide.json` to point to the described input files and change parameters
4. Run fitDNM snakemake 

### LSF submisison example
`bsub  -R 'rusage[mem=50GB]' -n 1 -a 'docker(docker/dockerfile)' /opt/conda/envs/snakemake/bin/snakemake -s /path/to/fitDNM/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1 `


# Input files:
- CADD files (v1.6) and md5sums:
  - `whole_genome_SNVs.tsv.gz  faaa80ef3948cf44e56a3629a90cdaaa`  
  - `whole_genome_SNVs.tsv.gz.tbi  4843cab24dd4992bb0cc5f1a7ebc807a`
- `mutation_rate_by_trinucleotide_matrix.txt ca2faad78f1055c266b5a8451bebf1cb`

### File created in pipeline:
- Input file: Bedfile containing the regions you want to investigate
- `.CADD.txt`: comprehensive list of CADD scores, this table is in a different format then the one downloaded from the CADD website, where instead of each possible change for a given nucleotide is an entry, each nucleotide is an entry and the changes are columns. The *PHRED* score is used for fitDNM, not the raw score.
- `.lis`: list of mutations in your region of interest, columns should be chr, pos, ref, alt, gene
- `mutation_rate_by_trinucleotide_matrix.txt`: The frequencies of mutations for each given trinucleotide. In this file what youre looking at is the frequency for *second* nucleotide  being changed
- `.mu.lis`: utilizes the trinucleotide mutation rate frequencies to calculate the mutation rate for every possible change


### dockerfile requirements:
- tabix
- Snakemake
- R and the following R packages:
  - foreach
  - iterators
  - doParallel


### Known issues:
  - When running on more than core the work flow errors out frequently in rule get_region_mutation or rule calculate_trinucleotide_mutations and returns an error saying the index is out of range.

### To do list:
- [] Update rule get_region_mutation to use tabix instead of python script to improve speed
  - [] zip and index the mutation file  
- [] Add in filter to seperate proband and sibling mutations
- [] fix mutation report file, right now it adds a header for every mutation file created
