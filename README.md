# Agalepto - MSc thesis script library

Collection of R and Bash scripts for analysis and visualization of reduced representation sequencing data (nextRAD). These scripts were used to assess the "true" diversity within the mesophotic coral genera *Leptoseris* and *Agaricia* for the MSc thesis:  Undescribed coral diversity at mesophotic depths: a phylogenomic assessment of *Agaricia* and *Leptoseris*. Although all scripts found here are functional, some of them are hardcoded and might need some small editing to be used for other data sets. This repository represents therefore a work in progress.

### admixture.sh 

Runs admixture on bash for RADseq data, and output to a `.csv` file containing the membership values of each group.

### agalepto_map_fig1.R

Creates maps of different geographic regions depending on the latitude/longitude given. These maps are created in R using R studio and can be exported to use with other softwares.

### denovo_DAPC.R

Runs denovo DAPC analysis (unsupervised, without assigning populations) for replicate `.vcf` files and output to a `.csv` file containing the membership values of each group. This analysis is created in R using R studio and can be exported to use with other softwares.

### denovo_DAPC_plot.R

Runs denovo DAPC analysis (unsupervised, without assigning populations) for replicate `.vcf` files and output to a `.csv` file containing the membership values of each group, and in addition it can create plots to visualize the membership values. This analysis is created in R using R studio and can be exported to use with other softwares.

### fig2-leptoseris_main_figure.R

Create a composition figure for *Leptoseris* dataset, including plotting a species tree, DAPC results, metadata, TSNE results. This figure is created in R studio and can be exported to use with other softwares.

### fig3_leptoseris.R

Create a composition figure for *Leptoseris* dataset, including plotting a species tree, DAPC results, metadata, BFD results. This figure is created in R studio and can be exported to use with other softwares.

### tanglegram_raxml_cox.R

Create tangle tree of RAD vs COX1 marker to compare two phylogenies using different genetic markers. This figure is created in R studio and can be exported to use with other softwares.

