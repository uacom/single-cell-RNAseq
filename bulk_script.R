#!/usr/bin/env Rscript

library(Rsubread)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")  # get reference data for mm10
library(BSgenome.Mmusculus.UCSC.mm10)

file_names = c("FAD1190", "FAD1193", "FAD1219", "wt1191", "wt1213", "wt1214")

#Data pre-processing: import and find all of the fastq.gz files
fastq1190.files = list.files(path = "./usftp21.novogene.com/raw_data/FAD1190", pattern = ".fq.gz$", full.names = TRUE)
fastq1193.files = list.files(path = "./usftp21.novogene.com/raw_data/FAD1193", pattern = ".fq.gz$", full.names = TRUE)
fastq1219.files = list.files(path = "./usftp21.novogene.com/raw_data/FAD1219", pattern = ".fq.gz$", full.names = TRUE)
fastq1191.files = list.files(path = "./usftp21.novogene.com/raw_data/wt1191", pattern = ".fq.gz$", full.names = TRUE)
fastq1213.files = list.files(path = "./usftp21.novogene.com/raw_data/wt1213", pattern = ".fq.gz$", full.names = TRUE)
fastq1214.files = list.files(path = "./usftp21.novogene.com/raw_data/wt1214", pattern = ".fq.gz$", full.names = TRUE)

fastq.files = c(fastq1190.files, fastq1193.files, fastq1219.files, fastq1213.files, fastq1214.files, fastq1191.files)

# write out reference sequence in FASTA format:
mm10 = BSgenome.Mmusculus.UCSC.mm10

Rsubread::buildindex(basename="mm10",
                      reference="GRCm38.p6.genome.fa")

