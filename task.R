# task.R
# project title: 'BiGCaT Coding Assignment'
# author: Woosub Shin
# final update: 16 November 2019
###############################################
# This scripts performs three tasks for BiGCaT.
# First task is to link a list of SNPs 
# to appropriate genes using biomart. 
# Second task is to create a SNP-gene network
# on Cytospace. 
# Third, and final task is to extend 
# the SNP-gene network by WikiPathway linkset, 
# using CyTargetLinker. 
#
# Each task section is marked and commented. 
###############################################


source("constants.R")
setwd(MAIN_DIR)

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
colnames(edges_df) <- c(SOURCE, TARGET)

# create the network
createNetworkFromDataFrames(nodes_df,
                            edges_df,
                            title=SNP_GENE_NETWORK,
                            directed=FALSE,
                            collection=GRAPHNEL_NETWORKS)

#########################################################
#   Task 3: Extend Network using Wikipathways linkset   #
#########################################################

# adapted from github:cytargetlinker-automation/R-automation/UseCase2/UseCase2.Rmd
wp <- file.path(getwd(), LINKSETS, WIKIPATHWAY_FILE)
CTLextend.cmd = paste('cytargetlinker extend idAttribute=Ensemble linkSetFiles="',
                      wp,
                      '" network="',
                      SNP_GENE_NETWORK,
                      '" direction=BOTH',
                      sep="")
commandsRun(CTLextend.cmd)
layoutNetwork()

# save and export file 
png_file <- file.path(getwd(), "snp_gene_network_extended.png")
exportImage(png_file, PNG)

# apply styles (again, adapted from above)
vizstyle_file <- file.path(getwd(), VIZSTYLES, VIZ_USECASE2)
LoadStyle.cmd = paste('vizmap load file file="', vizstyle_file,'"', sep="")
commandsRun(LoadStyle.cmd)
ApplyStyle.cmd = 'vizmap apply styles="CTL Gene Pathway Network"'
commandsRun(ApplyStyle.cmd)

# save and export file 
png_file <- file.path(getwd(), "snp_gene_network_extended2.png")
exportImage(png_file, PNG)
