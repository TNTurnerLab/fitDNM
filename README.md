# fitDNM
### Evin Padhi
### Last updated: 3/8/21
### Washington University in St. Louis
### School of Medicine, Department of Genetics
### Laboratory of Dr. Tychele Turner
### Snakefile for running fitDNM and all preprocessing steps
### Contact: evin.padhi@wustl.edu
### http://turnerlab.wustl.edu


### Purpose:
- To create a workflow that can take in a set of genomic coordinates, variants, CADD scores, and mutations rates to run through fitDNM in a high-throughput manner

### fitDNM:
- Originally developed by the Allen lab (http://people.duke.edu/~asallen/Software.html) in Jiang et al 2015, *Am. J. Hum. Genet.*  (https://www.cell.com/ajhg/fulltext/S0002-9297(15)00277-3)


### How to run:
1. Download CADD scores and check md5sums
2. Build dockerfile and push to dockerhub
3. Modify .json  `fitDNM_genome_wide.json` to point to the described input files and change parameters
4. Run fitDNM snakemake. Once finished running, two files should be generated, `.fitDNM.report` and `.muts.report` both detailed below

### LSF submission example
`bsub  -R 'rusage[mem=10GB]' -n 1 -a 'docker(docker/dockerfile)' /opt/conda/envs/snakemake/bin/snakemake -s /path/to/fitDNM/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1 `


### Local run example:
- Assuming that your data and fitDNM snakemake are in separate paths and that you have built the dockerfile and pushed it to dockerhub, you could run the following. Just be sure to update the paths in `fitDNM_genome_wide.json` and the path to the configfile in `fitDNM_genome_wide.smk` to reflect what you are mounting in the docker run command.
`docker run -v "/path/to/fitDNM/fitDNM_snakemake:/fitDNM_snakemake" -v "/path/to/data:/data" user/fitDNM_snakemake:latest /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_snakemake/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1`


# Input files:
- CADD files and md5sums https://cadd.gs.washington.edu/download :
  - `whole_genome_SNVs.tsv.gz  faaa80ef3948cf44e56a3629a90cdaaa`  
  - `whole_genome_SNVs.tsv.gz.tbi  4843cab24dd4992bb0cc5f1a7ebc807a`
     Note we only support CADD score release v1.6 for hg38, description: All possible SNVs of GRCh38/hg38
- `mutation_rate_by_trinucleotide_matrix.txt ca2faad78f1055c266b5a8451bebf1cb`
- Bed file that using standard format of `chr +'\t' + start + '\t' + annotation + '\n'` of the regions of interest
  - _each entry must have a unique annotation_
- Mutation file, to integrate automatically into workflow the mutation file must be ordered in the following
  `holder_column + " " + chromosome + " " +  position + " " + reference + " " + alternate + '\n'`


### Files created in pipeline:
- `.CADD.txt`: comprehensive list of CADD scores, this table is in a different format then the one downloaded from the CADD website, where instead of each possible change for a given nucleotide is an entry, each nucleotide is an entry and the changes are columns. The _PHRED_ score is used for fitDNM, not the raw score.
- `.lis`: list of mutations in your region of interest, columns should be chr, pos, ref, alt, gene
- `mutation_rate_by_trinucleotide_matrix.txt`: The frequencies of mutations for each given trinucleotide. In this file what you're looking at is the frequency for *second* nucleotide  being changed
- `.mu.lis`: utilizes the trinucleotide mutation rate frequencies to calculate the mutation rate for every possible change
- `.fitDNM.report` contains the results of fitDNM for all elements in the bedfile, which should consist of 8 columns for elements that have SNVs
- `.muts.report` summarizes the mutations in each element.


### Dockerfile requirements:
- Tabix (http://www.htslib.org/doc/tabix.html)
- Snakemake (https://snakemake.readthedocs.io/en/stable/)
- Bedtools (https://bedtools.readthedocs.io/en/latest/)
- R and the following R packages:
  - foreach (https://cran.r-project.org/web/packages/foreach/foreach.pdf)
  - iterators (https://cran.r-project.org/web/packages/iterators/iterators.pdf)
  - doParallel (https://cran.r-project.org/web/packages/doParallel/doParallel.pdf)


### Known issues:
- tabix continually runs on a small amount of files on runs that have a large amount of elements scanned (fixed 3/18/21)
