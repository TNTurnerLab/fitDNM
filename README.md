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


### File requirements (The snakemake should take care of everything but heres an overview):
- Input file: Bedfile containing the regions you want to investigate
- `.CADD.txt`: comprehensive list of CADD scores, this table is in a different format then the one downloaded from the CADD website, where instead of each possible change for a given nucleotide is an entry, each nucleotide is an entry and the changes are columns. The *PHRED* score is used for fitDNM, not the raw score.
- `.lis`: list of mutations in your region of interest, columns should be chr, pos, ref, alt, gene
- `mutation_rate_by_trinucleotide_matrix.txt`: The frequencies of mutations for each given trinucleotide. In this file what youre looking at is the frequency for *first* nucleotide  being changed
- `.mu.lis`: utilizes the trinucleotide mutation rate frequencies to calculate the mutation rate for every possible change


### software for dockerfile:
- tabix
- Snakemake
- R and the following R packages:
  - foreach
  - iterators
  - doParallel
