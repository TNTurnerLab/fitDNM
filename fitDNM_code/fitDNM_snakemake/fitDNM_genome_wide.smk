# Evin Padhi
# 1/30/21
# Washington University in St. Louis
# School of Medicine, Department of Genetics
# Laboratory of Dr. Tychele Turner
# Snakefile for running fitDNM and all preprocessing steps
# Contact: evin.padhi@wustl.edu


import os

configfile: "/fitDNM_code/fitDNM_snakemake/fitDNM_genome_wide.json"

# give absolute path to all files, input file must be a bed file
# with the columns being chr, start, stop, region name
REGIONS_OF_INTEREST = config['regions_of_interest']
CADD_SCORE_FILE = config['cadd_score_file']
MUTATION_CALLS = config['mutation_calls']
TRINUCLEOTIDE_MUT_RATE = config['trinucleotide_mut_rate']
MALES = config['males']
FEMALES = config['females']
FITDNM_PATH = config['fitDNM_R_path']
SADDLE_POINT_PATH = config['saddle_point_path']
R_SCRIPT_PATH = config['transform_cadd_scores_script_path']

# creates file names for output files
OUTPUT_PREFIX = REGIONS_OF_INTEREST.strip().split('/')[-1].split('.')[0]
OUTPUT_FILE_FITDNM = str(OUTPUT_PREFIX) + ".fitDNM.report"
OUTPUT_FILE_MUTS = str(OUTPUT_PREFIX) + ".muts.report"

BED_FILE_NAME = REGIONS_OF_INTEREST.split('/')[-1]

# creates a bed file for each unique entry in the regions list and also gets these names to fill wildcards
# by adding them to a list
region_annotations = []
for region in open(REGIONS_OF_INTEREST):
    region_processed = region.strip().split()
    region_name = region_processed[3]
    chrom_only = region_processed[0].split('r')[1]
    region_annotations.append(region_name)
    out_file_name = str(region_name) + ".tabix.bed"
    with open(out_file_name,"w") as out:
        out.write(chrom_only + '\t' + region_processed[1] + '\t' +region_processed[2])



rule all:
  input:
    OUTPUT_FILE_FITDNM,
    OUTPUT_FILE_MUTS,
    CADD_SCORE_FILE,
    MUTATION_CALLS,
    TRINUCLEOTIDE_MUT_RATE,
    REGIONS_OF_INTEREST,
    "intersected_%s" % BED_FILE_NAME



# Because the trinucleotide mutation rate calculates the mutation rate of the middle nucleotide in
# a trinucleotide we have to get a sequence that is plus one to each end because otherwise you lose
# the end bases in the final analysis. The best way to do this is to just grab the sequence from the CADD score
# file again.
rule get_sequence_plus_two_and_cadd_scores:
    input:region="{genomic_region}.tabix.bed", scores="%s" % CADD_SCORE_FILE
    params:script_path=R_SCRIPT_PATH,region="{genomic_region}"
    output: sequence_plus_two_final=temp("{genomic_region}_plus_two.fasta"),final_cadd=temp("{genomic_region}.final.CADD.txt"),final_fasta=temp("{genomic_region}.fasta")
    shell:"""
    cat {input.region} | tabix {input.scores} -R - | Rscript {params.script_path}/transform_CADD_scores.R {params.region} {output.final_cadd} {output.final_fasta}
    cat {input.region} |  awk -v s=0 '{{print $1, int($2)-1,int($3)+1}}' | tabix {input.scores} -R - | awk '!array[$2]++'  | awk '{{print $3}}' | tr -d "\n\r" > {output.sequence_plus_two_final}
    """

# Uses bedtools to number of variatns and variant position and status for each
# entry in the bedfile
rule intersect:
    input:regions=REGIONS_OF_INTEREST ,DNMs=MUTATION_CALLS
    output:intersected=temp("intersected_%s" % BED_FILE_NAME)
    shell:"""
    cat {input.DNMs} | awk '{{print $2,int($3),int($3)+1,$4,$5}}' OFS='\t' | /opt/conda/bin/bedtools  intersect -a {input.regions} -b - -wa -wb | /opt/conda/bin/bedtools groupby -i - -g 1-4 -c 6,8,9 -o collapse,collapse,collapse > {output.intersected}
    """

# Loops through intersected file and creates a file for each three files for each
# entry in the bed file, 1. a list of all mutations, 2. a list of all SNVs, and
# 3. a list of all MNVs
rule get_region_mutation:
    input:"intersected_%s" % BED_FILE_NAME
    output:temp(expand("{genomic_region}.lis",genomic_region=region_annotations)),temp(expand("{genomic_region}.snv.lis",genomic_region=region_annotations)),temp(expand("{genomic_region}.mnv.lis",genomic_region=region_annotations))
    run:
        import sys

        def create_list_files(line):
            # extracts needed information from line
            line_entry = line.strip().split('\t')
            ID,variant_position,ref,alt = line_entry[3],line_entry[4],line_entry[5],line_entry[6]
            out_prefix = str(ID)

            #gets the number of the chromosome
            number_only = line_entry[0].split('r')[1]


            # checks for number of variants and if more than one variant
            # creates lists to hold the information of each variant
            if len(variant_position.split(',')) > 1 :
                variant_coordinates = variant_position.split(',')
                ref_alleles = ref.split(',')
                alt_alleles = alt.split(',')
                with open(out_prefix+".lis",'w') as total_out,open(out_prefix+".snv.lis",'w') as snv_out,open(out_prefix+".mnv.lis",'w') as mnv_out:

                    total_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                    snv_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                    mnv_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file

                    for position_entry,ref_entry,alt_entry in zip(variant_coordinates,ref_alleles,alt_alleles):
                        total_out.write(number_only + '\t' + position_entry + '\t' + ref_entry + '\t' + alt_entry + '\t'  + ID + '\n')
                        if len(ref_entry) > 1 or len(alt_entry) > 1:
                            mnv_out.write(number_only + '\t' + position_entry + '\t' + ref_entry + '\t' + alt_entry + '\t'  + ID + '\n')
                        else:
                            snv_out.write(number_only + '\t' + position_entry + '\t' + ref_entry + '\t' + alt_entry + '\t'  + ID + '\n')

            # does same as above but with only variant present
            if len(variant_position.split(',')) == 1:

                with open(out_prefix+".lis",'w') as total_out,open(out_prefix+".snv.lis",'w') as snv_out,open(out_prefix+".mnv.lis",'w') as mnv_out:

                    total_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                    snv_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                    mnv_out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                    total_out.write(number_only + '\t' + variant_position + '\t' + ref + '\t' + alt + '\t'  + ID + '\n')

                    if len(ref) >1 or len(alt) >1:
                        mnv_out.write(number_only + '\t' + variant_position + '\t' + ref + '\t' + alt + '\t'  + ID + '\n')
                    else:
                        snv_out.write(number_only + '\t' + variant_position + '\t' + ref + '\t' + alt + '\t'  + ID + '\n')
                        statement = "No MNVs in " + str(ID)
                        mnv_out.write(statement)

        for line in open(input[0]):
            create_list_files(line)


rule calculate_trinucleotide_mutations:
    input: mut_rate=TRINUCLEOTIDE_MUT_RATE,sequence="{genomic_region}_plus_two.fasta",coordinates="{genomic_region}.tabix.bed",wait_file="{genomic_region}.lis",wait2_file="{genomic_region}.snv.lis",wait3_file="{genomic_region}.mnv.lis"
    params: annotation="{genomic_region}"
    output:temp("{genomic_region}.mu.lis")
    run:
        # Initializes dictionary to hold the mutation rates for different trinucleotides
        # where each trinucleotide is a key and the mutation rate of the middle nucleotide
        # are the values
        # order of values is: A T C G
        mutation_rates = {}

        #initializes counter to skip header
        counter = 0

        # opens mutation rate file and loops through it
        # and uses the counter to skip the header
        for trinucleotide_entry in open(input[0]):
            counter += 1
            if counter > 1:

                # gets trinucleotide and all corresponding mutation rates
                processed_entry = trinucleotide_entry.strip().split()
                mutation_rates[processed_entry[0]] = [processed_entry[1],processed_entry[2],processed_entry[3],processed_entry[4]]


        # gets coordinates from bed file
        coordinates = [[entry for entry in region.strip().split('\t')] for region in open(input[2])]

        # opens output mutation file and writes header
        # and initializes position counter
        with open(output[0],"w") as out:
            out.write("Annotation" + "\t" + "chr" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "A" + "\t" + "T" + "\t" + "C" + "\t" + "G" + '\n')
            position_counter = -1

            # loops through file to get fasta and
            # then loops through sequence, adding to
            # the sequence counter for each nucleotide
            for line in open(input[1]):
                for x in range(len(line)):
                        position_counter += 1

                        # calculates current position and extracts the trinucleotide
                        # and checks to make sure it is a trinucleotide
                        current_position = int(coordinates[0][1]) + position_counter
                        trinucleotide = line[x:x+3]
                        if len(trinucleotide) < 3:
                            break

                        # if the sequence is a trinucleotide writes to file
                        # with all the calculated rates
                        else:
                            if trinucleotide in mutation_rates:
                                out.write( str(params.annotation) + '\t'+ str(coordinates[0][0]) + '\t'+ str(current_position) + '\t'+ line[x+1] +  '\t' + mutation_rates[trinucleotide][0] + '\t' + mutation_rates[trinucleotide][1] + '\t' + mutation_rates[trinucleotide][2] + '\t' + mutation_rates[trinucleotide][3] + '\n')



# Runs fitDNM using all of the previously created files
rule run_fitDNM:
    input: mutation_list="{genomic_region}.snv.lis",all_muts="{genomic_region}.mu.lis",CADD_file="{genomic_region}.final.CADD.txt",tabix="{genomic_region}.tabix.bed"
    output:init=temp("{genomic_region}.fitDNM_CADD.txt")
    params: number_of_males=MALES,number_of_females=FEMALES, R_code_path="%s" % FITDNM_PATH, saddle_point_path=SADDLE_POINT_PATH, annotation="{genomic_region}"
    shell:"""
    Rscript {params.R_code_path}fitDNM.CADD.R {params.number_of_males} {params.number_of_females} {input.mutation_list} {input.all_muts} {input.CADD_file} {output.init} {params.saddle_point_path} {params.annotation}
    """

# Adds all fitDNM output into one file
rule clean:
    input: expand("{genomic_region}.fitDNM_CADD.txt",genomic_region=region_annotations)
    output:OUTPUT_FILE_FITDNM
    shell:"""
    echo "geneID nsample.adj nsnv nsnv.analysis ndenovo score pvalue.fitDNM pvalue.Poisson" > {output}
    cat {input} | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8}}' >> {output}
    """

# Combines all mutation files into one file
rule combine_mutations:
    input: mut_files=expand("{genomic_region}.lis",genomic_region=region_annotations),wait_files=OUTPUT_FILE_FITDNM
    output: OUTPUT_FILE_MUTS
    shell:"""
    echo "chr position Ref Alt annotation" > {output}
    cat {input.mut_files} | awk '{{print $1,$2,$3,$4,$5}}' >> {output}
    rm -rf *.tabix.bed
    """
