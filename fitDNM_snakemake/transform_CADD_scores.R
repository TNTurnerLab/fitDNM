# Evin Padhi
# 1/30/21
# Washington University in St. Louis
# School of Medicine, Department of Genetics
# Laboratory of Dr. Tychele Turner
# Snakefile for running fitDNM and all preprocessing steps
# Contact: evin.padhi@wustl.edu

library(tidyr)

args <-commandArgs(trailingOnly=TRUE)
sequence_file_name <- toString(args[3])
scores_file_name <- toString(args[2])
region <- toString(args[1])


in_file <- file("stdin")
CADD_file <- read.table(in_file, sep = "\t",header = FALSE)

#CADD_file <- read.table(snakemake@input[[1]], sep = "\t",header = FALSE)
CADD_file$V5 <- NULL
CADD_file_processed <- tidyr::pivot_wider(CADD_file,names_from= V4, values_from = V6,values_fill = 0)
colnames(CADD_file_processed)[1] <- "chr"
colnames(CADD_file_processed)[2] <- "pos"
colnames(CADD_file_processed)[3] <- "Ref"
#CADD_file_processed$Gene <- rep(snakemake@params[[2]],length(CADD_file_processed$chr))
CADD_file_processed$Gene <- rep(region,length(CADD_file_processed$chr))



rearranged_CADD_file <- CADD_file_processed[,c("Gene","chr","pos","Ref","A","T","C","G")]




sequence <- paste(rearranged_CADD_file$Ref, collapse = "")
write(sequence, file=sequence_file_name)

#write(sequence, file=snakemake@output[[2]])

#write.table(rearranged_CADD_file,file=snakemake@output[[1]], sep='\t',col.names=TRUE,quote=FALSE,row.names = FALSE)
write.table(rearranged_CADD_file,file=scores_file_name, sep='\t',col.names=TRUE,quote=FALSE,row.names = FALSE)
