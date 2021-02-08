# Evin Padhi
# 1/30/21
# Washington University in St. Louis
# School of Medicine, Department of Genetics
# Laboratory of Dr. Tychele Turner
# Snakefile for running fitDNM and all preprocessing steps
# Contact: evin.padhi@wustl.edu


import os

if config == {}:
  configfile: "/fitDNM/code/fitDNM_snakemake/fitDNM_genome_wide.json"

# give absolute path to all files, input file must be a bed file
# with the columns being chr, start, stop, region name
REGIONS_OF_INTEREST = config['regions_of_interest']
CADD_SCORE_FILE = config['cadd_score_file']
MUTATION_CALLS = config['mutation_calls']
TRINUCLEOTIDE_MUT_RATE = config['trinucleotide_mut_rate']
SAMPLE_SIZE = config['sample_size']
FITDNM_PATH = config['fitDNM_R_path']
SADDLE_POINT_PATH = config['saddle_point_path']

# creates file names for output files
OUTPUT_PREFIX = REGIONS_OF_INTEREST.strip().split('/')[-1].split('.')[0]
OUTPUT_FILE_FITDNM = str(OUTPUT_PREFIX) + ".fitDNM.report"
OUTPUT_FILE_MUTS = str(OUTPUT_PREFIX) + ".muts.report"


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
    "done.txt"



# condeses CADD scores to a format where each position is only one entry, instead of each change being an entry
# also gets sequence for mutation rate calculation
rule transform_cadd_scores:
    input:"{genomic_region}.initial.CADD.txt"
    group:"produce_files"
    output:temp("{genomic_region}.final.CADD.txt"),temp("{genomic_region}.fasta")
    params: "{genomic_region}"
    script: "transform_CADD_scores.R"


# Becasue the trinucleotide mutation rate calculates the mutation rate of the middle nucleotide in
# a trinucleotide we have to get a sequence that is plus one to each end because otherwise you lose
# the end bases in the final analysis. The best way to do this is to just grab the sequence from the CADD score
# file again.
rule get_sequence_plus_two:
    input:region="{genomic_region}.tabix.bed", scores="%s" % CADD_SCORE_FILE
    group:"produce_files"
    output: sequence_plus_two=temp("{genomic_region}_plus_two_init.fasta"), bed_plus_two=temp("{genomic_region}.plustwo.bed"), sequence_plus_two_final=temp("{genomic_region}_plus_two.fasta")
    shell:"""
    cat {input.region} |  awk -v s=0 '{{print $1, int($2-1),int($3+1)}}' > {output.bed_plus_two}
    cat {input.region} |  awk -v s=0 '{{print $1, int($2-1),int($3+1)}}' | tabix {input.scores} -R {output.bed_plus_two} | awk '!array[$2]++' > {output.sequence_plus_two}
    cat {output.sequence_plus_two} | awk '{{print $3}}' | tr -d "\n\r" > {output.sequence_plus_two_final}
    """


# uses tabix to get CADD scores
rule get_cadd_scores:
    input: region="{genomic_region}.tabix.bed", scores="%s" % CADD_SCORE_FILE
    group:"produce_files"
    output: initial=temp("{genomic_region}.initial.CADD.txt")
    shell: """
    tabix {input.scores} -R {input.region}  > {output.initial}
    """

# gets number of mutations in regions seperately
# mutation calls must be in a format where the chromosome is the second column
# the position is the third column, reference allele fourth
# and alternate allele fifth
# the first column can be anything but must be there
rule get_region_mutation:
    input:regions="{genomic_region}.tabix.bed" ,DNMs=MUTATION_CALLS
    params: annotation="{genomic_region}"
    group:"produce_files"
    output:temp("{genomic_region}.lis")
    run:
        import sys

        #gets coordinates from bed file, returns list of lists
        region = [[entry for entry in region.strip().split('\t')] for region in open(input[0])]

        # creates entry that contains chr# to look for matches in DNM file
        chromosome_entry = 'chr' + region[0][0]

        # opens output file and writes header
        # then loops through bed file and scans each
        # entry for mutations
        with open(output[0], "w") as out:
            out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
            for line in open(input[1]):
                line_p = line.strip().split('\t')

                # searches for mutations by first looking for entries
                # that match chromosomes then checkes in position is in range
                # and writes mutation to file if there is a match
                if line_p[1] == chromosome_entry: #
                    if int(line_p[2]) in range(int(region[0][1]),int(region[0][2])):
                        number_only = line_p[1].split('r')[1]
                        out.write(number_only + '\t' + line_p[2] + '\t' + line_p[3] + '\t' + line_p[4] + '\t'  + str(params.annotation) + '\n')


# creates two files from the mutation list,
# a file that holds SNVs andn a file that holds MNVs
rule filter_multinucleotide_variants:
    input:"{genomic_region}.lis"
    params:annotation="{genomic_region}"
    output:filtered=temp("{genomic_region}.snv.lis"),mnv_list=temp("{genomic_region}.mnv.lis")
    run:
        # initializes lists to hold SNVs and MNVs
        # that will later be written to files
        SNVs = []
        MNVs = []

        # counts number of MNVs
        # if 0 it will trigger a print statement
        # verifying no MNVs in region
        MNV_counter =0

        #initializes line_counter to skip first header line
        line_counter = 0

        #loops through mutation file
        for entry in open(input[0]):
            line_counter += 1
            if line_counter > 1:
                line_p = entry.strip().split('\t') # processes entries
                if len(line_p[2]) > 1 or len(line_p[3]) > 1: # checks lenght of reference and alternate allele and if longer than 1, adds to MNV list
                    MNVs.append(entry) #adds entry to MNV list
                    MNV_counter += 1
                else:
                    SNVs.append(entry) # appends mutation to SNV list if both ref and alt alleles are equal to 1

        #writes SNV file
        with open(output[0],"w") as out:
            out.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
            for entry in SNVs: # loops through SNVs list
                out.write(entry) # writes SNV to file

        #writes MNV file
        with open(output[1],"w") as out2:
            if MNV_counter < 1:
                statement = "No MNVs in " + str(params.annotation)
                out2.write(statement)
            else: # if MNVs do exist
                out2.write('chr' + '\t' + 'position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' +  'annotation' + '\n') #adds header to file
                for entry in MNVs: #loops through MNV list
                    out2.write(entry) # writes MNV  to file


rule calculate_trinucleotide_mutations:
    input: mut_rate=TRINUCLEOTIDE_MUT_RATE,sequence="{genomic_region}_plus_two.fasta",coordinates="{genomic_region}.tabix.bed"
    params: annotation="{genomic_region}"
    group:"produce_files"
    output:temp("{genomic_region}.mu.lis")
    run:
        # initialize dictionary to hold the mutation rates for different trinucleotides
        # where each trinucleotide is a key and the mutation of rate of the middle nucletotide
        # to change are the values
        # order of values is: A T C G
        mutation_rates = {}

        #initializes counter to skip header
        counter = 0

        # opens mutation rate file and loops through it
        # and uses a counter to skip the header
        for trinucleotide_entry in open(input[0]):
            counter += 1
            if counter > 1:

                # gets mutation trinucleotide and all corresponding mutation rates
                processed_entry = trinucleotide_entry.strip().split()
                mutation_rates[processed_entry[0]] = [processed_entry[1],processed_entry[2],processed_entry[3],processed_entry[4]]


        #gets coordinates from bed file, returns list of lists, change that later
        coordinates = [[entry for entry in region.strip().split('\t')] for region in open(input[2])]

        # opens output mutation file and writes header
        # and initializes position counter
        with open(output[0],"w") as out:
            out.write("Annotation" + "\t" + "chr" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "A" + "\t" + "T" + "\t" + "C" + "\t" + "G" + '\n')
            position_counter = 1

            # loops through file to get fasta and
            # then loops through sequence, adding to
            # the sequence counter for nucleotide
            for line in open(input[1]):
                for x in range(len(line)):
                        position_counter += 1

                        # calculates current position and extracts the tri nucleotide
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
    input: mutation_list="{genomic_region}.snv.lis",all_muts="{genomic_region}.mu.lis",CADD_file="{genomic_region}.final.CADD.txt"
    output:temp("{genomic_region}.fitDNM_CADD.txt")
    params: sample_size=SAMPLE_SIZE, R_code_path="%s" % FITDNM_PATH, saddle_point_path=SADDLE_POINT_PATH, annotation="{genomic_region}"
    shell:"""
    Rscript {params.R_code_path}fitDNM.CADD.R {params.sample_size} 0 {input.mutation_list} {input.all_muts} {input.CADD_file} {output} {params.saddle_point_path} {params.annotation}
    """

# Adds all fitDNM output into one file
rule clean:
    input: expand("{genomic_region}.fitDNM_CADD.txt",genomic_region=region_annotations)
    output:OUTPUT_FILE_FITDNM
    shell:"""
    echo "geneID nsample.adj nsnv nsnv.analysis ndenovo score pvalue.fitDNM pvalue.Poisson" > {output}
    cat {input} | awk 'NR % 2 == 0 {{print $1,$2,$3,$4,$5,$6,$7,$8}}' >> {output}
    """

# Combines all mutations into one file
rule combine_mutations:
    input: expand("{genomic_region}.lis",genomic_region=region_annotations)
    output: OUTPUT_FILE_MUTS
    shell:"""
    echo "chr position Ref Alt annotation" > {output}
    cat {input} | awk ' NR>1 {{print $1,$2,$3,$4,$5}}' >> {output}
    """

#echo "chr position Ref Alt annotation" > {output}

rule delete_tabix:
    input:output1=OUTPUT_FILE_MUTS,output2=OUTPUT_FILE_FITDNM,tabix=expand("{genomic_region}.tabix.bed",genomic_region=region_annotations)
    output: "done.txt"
    shell:"""
    rm -rf {input.tabix}
    touch {output}
    """
