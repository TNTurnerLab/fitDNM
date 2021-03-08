
export LSF_DOCKER_VOLUMES="/storage1/fs1/tychele/Active/projects/fitDNM:/fitDNM"


while getopts "f:" opt; do
  case ${opt} in
    f) bedfile=$OPTARG;;
    esac
done

FILENAME=${bedfile##*/}
BASENAME=${FILENAME%%.*}

if [ -z "$BASENAME" ]
then
  echo "EXIT: please provide input bedfile with -f"
  exit 1
else
  mkdir -p "$BASENAME.fitDNM_run"
  cd "$BASENAME.fitDNM_run"
  cp $bedfile .
  bsub -q tychele -R 'span[hosts=1] rusage[mem=200GB]' -n 10 -a 'docker(evinpadhi/smk_r_tabix:V2.0_bedtools)' /opt/conda/envs/snakemake/bin/snakemake -s /fitDNM/code/fitDNM_snakemake/fitDNM_genome_wide.smk \
  -C regions_of_interest=$FILENAME \
  --cores 10 \
  --restart-times 4
fi
