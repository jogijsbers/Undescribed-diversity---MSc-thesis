# Run denovo DAPC for replicate VCF files and output to CSV
# Johanna Gijsbers, Feb 2021
# Copyright: Copyright (C) 2021 J. Gijsbers Alejandre
# License: GPL


## Dependencies ==============================================================
library("adegenet")
library("vcfR")
library("ggplot2")
library("dplyr")
library("svglite")

## Main code =================================================================
# Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
arg_pattern <- "agalepto_lepto_singlesnps*.vcf" # args [1]
arg_cluster_start <- 4 # args [2]
arg_cluster_end <- 12 # args [3]
no_clusters = 1 + arg_cluster_end - arg_cluster_start
# Obtain list of vcf files to summarise
vcf_files <- list.files(path=".", all.files = TRUE, no.. = TRUE,
                        full.names=FALSE, recursive=FALSE, include.dirs = FALSE)

vcf_files <- list.files(path=".", pattern = arg_pattern, 
                        full.names=TRUE, recursive=FALSE)

# Loop over files and run find.clusters for each K
assignment_groups <- list()
for(filename in vcf_files){
  vcf_genlight <- vcfR2genlight(read.vcfR(filename))
  filebase <- sub('\\..*$', '', basename(filename))
  replicate = substr(filebase, 24, 24) # specific to agalepto_aga_singlesnp
 # Run find.clusters for each k
  for(i in arg_cluster_start:arg_cluster_end){
    assignment = find.clusters(vcf_genlight, max.n.clust = 20, n.clust = i, 
                               n.pca = 120, stat = c("BIC", "AIC", "WSS"), n.iter = 100000, 
                               n.start = 1000, set.seed(2))
    run_name = paste(replicate, "_k", i, sep="")
    assignment_groups[[run_name]] = assignment$grp
  }
}

assignment_groups_df <- as.data.frame(assignment_groups)
DAPC_posterior <- list()
DAPC_indcoord <- list()

# Run DAPC using the groups assignment from find.clusters 
for(i in colnames(assignment_groups_df)){
  DAPC = dapc(vcf_genlight, assignment_groups_df[[i]], n.pca = 40, # 18 for Agaricia, 3 for undaria, 40 for Leptoseris N/3
                n.da = 11) # 10 for Agaricia, 3 for Undaria, 10 for Leptoseris
  run_name = paste(replicate, "k", i, sep="")
  DAPC_posterior[[run_name]] = DAPC$posterior
  DAPC_indcoord[[run_name]] = DAPC$ind.coord
  
  }

# Combine individual runs into single data frame and save as csv
output_indcoord_filename_combined <- as.data.frame(DAPC_indcoord)
output_posterior_filename_combined <- as.data.frame(DAPC_posterior)
write.csv(output_indcoord_filename_combined, file = 'agalepto_lepto_denovo_DAPC_indcoord_5c2.csv')
write.csv(output_posterior_filename_combined, file = 'agalepto_lepto_denovo_DAPC_posterior_5c2.csv')
