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
                   header=FALSE,
                   stringsAsFactors=FALSE,
                   col.names=c(SNP),
                   sep="\t")

# use snp Ensembl database
snpmart <- useEnsembl(biomart=SNP,
                      dataset=HOMO_SAPIENS_SNP)

# if want to see vailable attributes & functions
listAttributes(snpmart)
listFilters(snpmart)

# get associated genes for each snp
snp_to_gene <- getBM(attributes=c(SNP_ID,
                                  ENSEMBLE_GENE_ID),
                     filters=SNP_FILTER,
                     values=snps,
                     mart=snpmart)

#########################################################
#    Task 2: Create a SNP-Gene Network in Cytoscape     #
#########################################################

# test connection to Cytoscape
cytoscapePing ()
cytoscapeVersionInfo ()

# get a unique collection of genes
unique_genes <- unique(unlist(snp_to_gene$ensembl_gene_stable_id,
                              use.names=FALSE))

# create dataframe for nodes
id <- c(unique_genes,
        unlist(snps, use.names=FALSE))
Ensemble <- c(unique_genes,
              rep(NA, nrow(snp_to_gene)))
nodes_df <- data.frame(id, Ensemble)

# create dataframe for edges
edges_df <- snp_to_gene
colnames(edges_df) <- c("source", "target")

# create the network
createNetworkFromDataFrames(nodes_df,
                            edges_df,
                            title='SNP-Gene Network',
                            directed=FALSE,
                            collection='GraphNEL Networks')

# save and export file (optional; for presentation)
png_file <- file.path(MAIN_DIR, "snp_gene_network.png")
exportImage(png_file,'PNG')

#########################################################
#   Task 3: Extend Network using Wikipathways linkset   #
#########################################################





