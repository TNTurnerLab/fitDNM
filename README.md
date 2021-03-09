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



- CADD files and their corresponding md5sums https://cadd.gs.washington.edu/download :
  - `whole_genome_SNVs.tsv.gz  faaa80ef3948cf44e56a3629a90cdaaa` corresponds to `cadd_score_file` in config file  
  - `whole_genome_SNVs.tsv.gz.tbi  4843cab24dd4992bb0cc5f1a7ebc807a`
  -  Note we only support CADD score release v1.6 for hg38 (file description: All possible SNVs of GRCh38/hg38, US link) and ensure both the `.tbi` file and `.tsv.gz` are downlaoded in the same directory
- `mutation_rate_by_trinucleotide_matrix.txt ca2faad78f1055c266b5a8451bebf1cb` corresponds to `trinucleotide_mut_rate` in config file  
- Bed file that corresponds to the genomic regions to be tested in fitDNM
  - please use standard format of `chr +'\t' + start + '\t' + annotation + '\n'`, make sure there are no extra columns after the fourth column   
  - corresponds to `regions_of_interest` in config file
  - _each entry must have a unique annotation_
- Variant file, corresponds to `mutation_calls` in config file
  - to integrate automatically into workflow the mutation file must be ordered in the following
  `holder_column + " " + chromosome + " " +  position + " " + reference + " " + alternate + '\n'`
- Ensure that after downloading all files that all md5sums match those provided above

### Input files:
| File name | Source | MD5Sum | annotation in configfile | 
|-----------| -------|------- | -------------------------|
| whole_genome_SNVs.tsv.gz|  https://cadd.gs.washington.edu/download | faaa80ef3948cf44e56a3629a90cdaaa` | `cadd_score_file`| 
|whole_genome_SNVs.tsv.gz.tbi| https://cadd.gs.washington.edu/download |  4843cab24dd4992bb0cc5f1a7ebc807a | NA |
| mutation_rate_by_trinucleotide_matrix.txt | Here | ca2faad78f1055c266b5a8451bebf1cb | `trinucleotide_mut_rate` | 
| Variant file | User provided | NA| `mutation_calls` | 
| Bed file | User provided | NA | `regions_of_interest` | 

All CADD score files can be downloaded from  https://cadd.gs.washington.edu/download using All possible SNVs of GRCh38/hg38 US link. Currently we oinly support release v1.6 for GRCh38/hg38  or alternatively use the wget commands below and be sure to check MD5sums after downloading to ensure the download was sucessful 
```
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

#### Formatting for user provided inputs:
For the user provided bed file please format using standard format with an annotation column for each entry:  
`chr +'\t' + start + '\t' + annotation + '\n'` <br/>  and make sure that each entry has a unique annotaiton 


### Files created in pipeline:
### Temporary files (fitDNM input)
 - `.CADD.txt`: comprehensive list of CADD scores, this table is in a different format then the one downloaded from the CADD website, where instead of each possible change for a given nucleotide is an entry, each nucleotide is an entry and the changes are columns. The _PHRED_ score is used for fitDNM, not the raw score.
 - `.lis`: list of mutations in your region of interest, columns should be chr, pos, ref, alt, gene
 - `mutation_rate_by_trinucleotide_matrix.txt`: The frequencies of mutations for each given trinucleotide. In this file what you're looking at is the frequency for *second* nucleotide  being changed
 - `.mu.lis`: utilizes the trinucleotide mutation rate frequencies to calculate the mutation rate for every possible change
### Final output
 - `.fitDNM.report` contains the results of fitDNM for all elements in the bedfile, which should consist of 8 columns for elements that have SNVs
 - `.muts.report` summarizes the mutations in each element.


## How to run:
1. Download CADD scores and check md5sums
2. Build dockerfile and push to dockerhub
3. Modify `fitDNM_genome_wide.json` and  to point to the described input files and change parameters
4. Modify `fitDNM_genome_wide.smk` to point to `fitDNM_genome_wide.json` path that is inside the docker image
5. Run fitDNM snakemake. Once finished running, two files should be generated, `.fitDNM.report` and `.muts.report` both detailed below


### LSF submission example
`export LSF_DOCKER_VOLUMES="/path/to/fitDNM_directory:/fitDNM"`
`bsub  -R 'rusage[mem=10GB]' -n 1 -a 'docker(docker/dockerfile)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1 `
 

### Local run example:
- Assuming that your data and fitDNM snakemake are in separate paths and that you have built the dockerfile and pushed it to dockerhub, you could run the following. Just be sure to update the paths in `fitDNM_genome_wide.json` and the path to the configfile in `fitDNM_genome_wide.smk` to reflect what you are mounting in the docker run command.
`docker run -v "/path/to/fitDNM/fitDNM_snakemake:/fitDNM_snakemake" -v "/path/to/data:/data" user/fitDNM_snakemake:latest /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1`

We provide a further example using a wrapper script to run the pipeline in `run_fitDNM_snake.sh` for an LSF system. For one to run this please change the `export LSF_DOCKER_VOLUMES` command to reflect where fitDNM is on your cluster and the respective bsub command to the memory and cpu requirements wanted and to pull the correct docker image. This command can be executed by running the following:
  `bash run_fitDNM_snake.sh -f /path/to/bedfile`

The script will then create a directory based of the name of your bed file and run fitDNM within the directory and generate the corresponding `.report` files

### Docker:
For those not familiar with docker please see https://docs.docker.com/get-started/overview/
- Example docker build command, first create a folder called `fitDNM_snakemake` that contains the Dockerfile from this repository in it. Then run
  1. `Docker build fitDNM_snakemake`
  2. `Docker images` to get the image ID
  3. `Docker tag <image_ID> <user_name>/fitDNM_snakemake:initial`
  4. `Docker push <user_name>/fitDNM_snakemake:initial`

#### Dockerfile requirements:
- Tabix (http://www.htslib.org/doc/tabix.html)
- Snakemake (https://snakemake.readthedocs.io/en/stable/)
- Bedtools (https://bedtools.readthedocs.io/en/latest/)
- R and the following R packages:
  - foreach (https://cran.r-project.org/web/packages/foreach/foreach.pdf)
  - iterators (https://cran.r-project.org/web/packages/iterators/iterators.pdf)
  - doParallel (https://cran.r-project.org/web/packages/doParallel/doParallel.pdf)


### Known issues:
- tabix continually runs on a small amount of files on runs that have a large amount of elements scanned (fixed 3/8/21)
