# fitDNM
### Evin Padhi
### Washington University in St. Louis
### Turner Lab

fitDNM was originally developed by the Allen lab (http://people.duke.edu/~asallen/Software.html) in Jiang et al 2015, *Am. J. Hum. Genet.*  (https://www.cell.com/ajhg/fulltext/S0002-9297(15)00277-3) to incoporate functional information in test of excess de novo mutational load. Here we've adapted the pipeline to utilize CADD scores instead of PolyPhen-2 scores in order to run in noncoding regions of the genome and implemented a scalable verision of the pipeline to test many elements at once. Given a bedfile that contains the regions of interest one wants to test for a significant excess of de novo mutations and the corresponding variants to use, this pipeline will output two summary files that contain the p values and scores calculated by fitDNM for each element in the bed file in the `.fitDNM.report` file and a summary of all mutations found in these genomic regions in the `.mutation.report` file

## Input files and download links:
| File name | Source | MD5Sum | annotation in configfile |
|-----------| -------|------- | -------------------------|
| `whole_genome_SNVs.tsv.gz` |  https://cadd.gs.washington.edu/download | faaa80ef3948cf44e56a3629a90cdaaa | `cadd_score_file`|
|`whole_genome_SNVs.tsv.gz.tbi` | https://cadd.gs.washington.edu/download |  4843cab24dd4992bb0cc5f1a7ebc807a | NA |
| `mutation_rate_by_trinucleotide_matrix.txt` | Here | ca2faad78f1055c266b5a8451bebf1cb | `trinucleotide_mut_rate` |
| Variant file | User provided | NA| `mutation_calls` |
| Bed file | User provided | NA | `regions_of_interest` |

All CADD score files can be downloaded from  https://cadd.gs.washington.edu/download using the All possible SNVs of GRCh38/hg38 US link, make sure to download both the score file and tabix index file or alternatively use the wget commands below and be sure to check MD5sums after downloading to ensure the download was sucessful.
```
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```
Currently we only support release v1.6 for GRCh38/hg38

__Formatting the user provided inputs:__<br>
Bed file: For the user provided bed file please format using standard bed file format with an extra annotation column for each entry:  
`chr +'\t' + start + '\t' + stop + '\t' + annotation + '\n'` __and make sure that each entry has a unique annotation__

Variant file: for the variant file to integrate automatically into the workflow the file must be ordered in the following way:
`holder_column + " " + chromosome + " " +  position + " " + reference + " " + alternate + '\n'` where the holder column can take on any value. They must also be on hg38

## Requirements and usage:

`Tabix`  (http://www.htslib.org/doc/tabix.html) <br>
`Snakemake`  (https://snakemake.readthedocs.io/en/stable/)<br>
`Bedtools` version 2.27.0 or above (https://bedtools.readthedocs.io/en/latest/)<br>
R and the following R packages:<br>
`foreach` (https://cran.r-project.org/web/packages/foreach/foreach.pdf)<br>
`iterators` (https://cran.r-project.org/web/packages/iterators/iterators.pdf)<br>
`doParallel` (https://cran.r-project.org/web/packages/doParallel/doParallel.pdf)<br>

__Outline of how to run:__ <br>
1. Clone repository using `git clone https://github.com/TNTurnerLab/fitDNM.git`
1. Download CADD scores and mutation rate file and check md5sums
2. Build dockerfile and push to dockerhub, alternatively pull from (insert link)
3. Format bed and variant file
4. Modify `fitDNM_genome_wide.json` to point to the described input files and change parameters also outlined below
5. Run fitDNM snakemake. Once finished running, two files should be generated, `.fitDNM.report` and `.muts.report`


To run change the paths in config file outlined below to reflect the actual paths of where the data is located (see section below for example setup)
```
{
  "mutation_calls": "/path/to/variants",
  "cadd_score_file": "/path/to/CADD_files/whole_genome_SNVs.tsv.gz",
  "trinucleotide_mut_rate": "/path/to/mutation_rate_by_trinucleotide_matrix.txt",
  "fitDNM_R_path": "/path/to/fitDNM_R_code/",
  "saddle_point_path": "/path/to/fitDNM_R_code/double_saddle_point_approx_8_7_2014.R",
  "males": "number_of_males_in_study",
  "females": "number_of_females_in_study",
  "transform_cadd_scores_script_path":"/path/to/fitDNM_snakemake",
  "regions_of_interest": "/path/to/bedfile"
}
```
To make running this pipeline easier when using docker, we suggest that you keep all input files in one parent directory that way you can run a single export or `-v` command for all of your input data  (e.g., an input_data directory that has different folders for CADD scores, variants, and regions of interest)


Note: We recommend creating a working directory when running this pipeline as it generates many different temporary files need to run fitDNM

__Running locally:__
If running locally with docker see the following example code and be sure to update the paths being mounted with `-v` to reflect the actual paths to the data and location of the fitDNM_snakemake code and the paths within the config file to reflect the mounted directories
```
docker run -v "/path/to/fitDNM/fitDNM_snakemake:/fitDNM_snakemake" -v "/path/to/data:/data" tnturnerlab/fitdnm_snakemake:V1.0 /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```

__Running on LSF:__
If running on an LSF server see below, be sure to follow the same steps and update all paths

```
export LSF_DOCKER_VOLUMES="/path/to/fitDNM:/fitDNM /path/to/input_data:/input_data"
bsub  -R 'rusage[mem=10GB]' -n 1 -a 'docker(tnturnerlab/fitdnm_snakemake:V1.0)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```
__Wrapper script:__ we also provide an example of a wrapper script that could be used on an LSF sever after updating the config file to point to all files needed except for the bedfile. The idea of this script is that it makes it easier to analyze different bed files  using the same set of variants by using `-f` instead of having to change the config file each time. To use this script first change the memory and cpu usage to the desired setting and update the paths and then run:
```bash run_fitDNM.sh -f /path/to/bed/file```
it will then create a directory for the run execution of the job and when finished contain all necessary output plus a copy of the bed file

### Example setups
Assuming all of our input files are in a directory called `input_data` and all of the code is in a directory called `fitDNM` where there is two sub-directories for the fitDNM stats code called `fitDNM_R_code` and a second directory for the snakemake and all preprocessing code called `fitDNM_snakemake` (which would be the setup upon pulling this repository), the config file and code to run on an LSF server would look like:
```
{
  "mutation_calls": "/fitDNM/code/input_data/input_data/variants_in.txt",
  "cadd_score_file": "/fitDNM/code/input_data/CADD_files/whole_genome_SNVs.tsv.gz",
  "trinucleotide_mut_rate": "/input_data/mutation_rate_by_trinucleotide_matrix.txt",
  "fitDNM_R_path": "/fitDNM/fitDNM_R_code/",
  "saddle_point_path": "/fitDNM/fitDNM_R_code/double_saddle_point_approx_8_7_2014.R",
  "males": "2666",
  "females": "0",
  "transform_cadd_scores_script_path":"/fitDNM/fitDNM_snakemake",
  "regions_of_interest": "/fitDNM/code/input_data/input_regions.bed"
}

export LSF_DOCKER_VOLUMES="/home/user/fitDNM_code:/fitDNM_code /home/user/input_data:/input_data"
bsub  -R 'rusage[mem=10GB]' -n 1 -a 'docker(docker/dockerfile)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_code/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```
Alternatively, running it locally and assuming the same file structure the command would look like:
```
docker run -v "/home/user/fitDNM_code:/fitDNM_code" -v "/home/user/data:/data" tnturnerlab/fitdnm_snakemake:V1.0 /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_code/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```

## Pipeline overview:

For each entry in the bedfile, this pipeline creates the following temporary files that are eventually used by fitDNM
 1. `<annotation>.CADD.txt`: comprehensive list of CADD scores, this table is in a different format then the one downloaded from the CADD website, where instead of each possible change for a given nucleotide is an entry, each nucleotide is an entry and the changes are columns. The _PHRED_ score is used for fitDNM, not the raw score.
 2. `<annotation>.lis`: list of mutations in your region of interest, columns should be chr, pos, ref, alt, gene
 3. `<annotation>.mu.lis`: utilizes the trinucleotide mutation rate frequencies to calculate the mutation rate for every possible change

After generating all neccesarry files and running fitDNM, it combines the fitDNM output of all entries into one file and all of the identified mutations into one file
 1. `.fitDNM.report` contains the results of fitDNM for all elements in the bedfile, which should consist of 8 columns for elements that have SNVs
 2. `.muts.report` summarizes the mutations in each element.



## Docker:
For those not familiar with docker please see https://docs.docker.com/get-started/overview/
- Example docker build command, first create a folder called `fitDNM_snakemake` that contains the Dockerfile from this repository in it. Then run
  1. `Docker build fitDNM_snakemake`
  2. `Docker images` to get the image ID
  3. `Docker tag <image_ID> <user_name>/fitDNM_snakemake:initial`
  4. `Docker push <user_name>/fitDNM_snakemake:initial`
