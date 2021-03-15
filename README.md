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
mkdir -p CADD_scores
cd CADD_scores
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



### Example setups
__Running an example analysis on the hs737 enhancer from Padhi et al 2020__ (https://www.biorxiv.org/content/10.1101/2020.08.28.270751v1):
After cloning the repository from github using `git clone https://github.com/TNTurnerLab/fitDNM.git` and creating the CADD_scores directory and downloading all files, please create a directory within `input_data` called `variants` to hold your variant call file. The directory set up should then resemble the one below, where the `.` represents the parent directory `fitDNM`

```
.
├── Dockerfile
├── README.md
├── fitDNM_code
│   ├── fitDNM_R_code
│   │   ├── double_saddle_point_approx_8_7_2014.R
│   │   └── fitDNM.CADD.R
│   └── fitDNM_snakemake
│       ├── fitDNM_genome_wide.json
│       ├── fitDNM_genome_wide.smk
│       └── transform_CADD_scores.R
├── input_data
│   ├── CADD_scores
│   │   ├── whole_genome_SNVs.tsv.gz
│   │   └── whole_genome_SNVs.tsv.gz.tbi
│   ├── input_regions
│   │   └── hs737.bed
│   ├── mutation_rate_by_trinucleotide_matrix.txt
│   └── variants
│       └── variants.txt
└── run_fitDNM_snake.sh
```


Using this setup, change the config file to the same as below:

```
{
  "mutation_calls": "/fitDNM/input_data/variants/variants.txt",
  "cadd_score_file": "/fitDNM/input_data/CADD_scores/whole_genome_SNVs.tsv.gz",
  "trinucleotide_mut_rate": "/fitDNM/input_data/mutation_rate_by_trinucleotide_matrix.txt",
  "fitDNM_R_path": "/fitDNM/fitDNM_R_code/",
  "saddle_point_path": "/fitDNM/fitDNM_code/fitDNM_R_code/double_saddle_point_approx_8_7_2014.R",
  "males": "2666",
  "females": "0",
  "transform_cadd_scores_script_path":"/fitDNM/fitDNM_code/fitDNM_snakemake",
  "regions_of_interest": "/fitDNM/input_data/hs737.bed"
}
```
Then to execute the code run one of the following

If running on an LSF server use the following:
```
export LSF_DOCKER_VOLUMES="/home/user/fitDNM_code:/fitDNM_code /home/user/input_data:/input_data"
bsub  -R 'rusage[mem=10GB]' -n 1 -a 'docker(tnturnerlab/fitdnm_snakemake:V1.0)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_code/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```

Alternatively, running it locally and assuming the same file structure the command would look like:

```
docker run -v "/home/user/fitDNM_code:/fitDNM_code" -v "/home/user/data:/data" tnturnerlab/fitdnm_snakemake:V1.0 /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM_code/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1
```

For fitDNM to run it is essential the correct paths are mounted in docker image. To determine the correct path please run `pwd` in the top level `fitDNM` directory to get the path. Then in both the `/home/user/fitDNM_code:/fitDNM_code` and `/home/user/input_data:/input_data` part of the code, replace the the lefthand side of the semicolon with the output from `pwd`. This should be done regardless of wether or not you are using an LSF server.


__Running on user files__
To run on user files, please place your variant file in the `variants` directory and your input bedfile into the `input_regions` directory and change the config file to point to these files.

__Wrapper script:__ we also provide an example of a wrapper script that could be used on an LSF sever after updating the config file to point to all files needed except for the bedfile. The idea of this script is that it makes it easier to analyze different bed files  using the same set of variants by using `-f` instead of having to change the config file each time. To use this script first change the memory and cpu usage to the desired setting and update the paths and then run:
`bash run_fitDNM.sh -f /path/to/bed/file`
it will then create a directory for the run execution of the job and when finished contain all necessary output plus a copy of the bed file


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
