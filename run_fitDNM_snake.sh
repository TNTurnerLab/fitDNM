export LSF_DOCKER_VOLUMES="/storage1/fs1/tychele/Active/projects/fitDNM:/fitDNM"


bsub -q tychele -R 'rusage[mem=50GB]' -n 1 -a 'docker(evinpadhi/smk_r_tabix:v1.0)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM/code/fitDNM_snakemake/fitDNM_genome_wide.smk --cores 1 
