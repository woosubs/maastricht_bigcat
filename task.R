# title: 'BiGCaT Interview Task'
# author: Woosub Shin
# date: 15 November 2019

setwd("~/Desktop/Study/Maastricht")
source("constants.R")

# load package and check
library(biomaRt)
library(reshape2)
library(RCy3)
listMarts()

#########################################################
#           Task 1: Link a Gene to Each SNP             #
#########################################################

# import names of SNPs
snps <- read.delim(SNP_FILE,
                   header = FALSE,
                   stringsAsFactors = FALSE,
                   col.names = c(SNP),
                   sep = "\t")

# use snp Ensembl database
snpmart <- useEnsembl(biomart = SNP,
                      dataset = HOMO_SAPIENS_SNP)

# if want to see vailable attributes & functions
listAttributes(snpmart)
listFilters(snpmart)

# get associated genes for each snp
snp_to_gene <- getBM(attributes = c(SNP_ID,
                                    ENSEMBLE_GENE_ID),
                     filters = SNP_FILTER,
                     values = snps,
                     mart = snpmart)


#########################################################
#    Task 2: Create a SNP-Gene Network in Cytoscape     #
#########################################################


