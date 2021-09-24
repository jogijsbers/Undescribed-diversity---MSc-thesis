## Author: JC. Gijsbers Alejandre
## Date: March 2021

### Create a composition figure for Leptoseris dataset, 
### including plotting a species tree, DAPC results, metadata, TSNE results 


## Dependencies ==============================================================
library("tidyverse")
library("ggplot2")
library("ggtree")
library("ggtreeExtra")
library("ggstance")
library("reshape2")
library("dplyr")
library("grid")
library("gtable")
library("phytools")
library("svglite")
library("aplot")
library("ggnewscale")
library("phylotools")
library("stringr")
library("cowplot")
library("ape")
library("ggpubr")
library("readr")

# Tree set-up
raxml <- read.tree("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4a - Maximum likelihood inference/agalepto_lepto_min10_raxml_suppTBE_4a.tree")
tips <- raxml$tip.label
dat <- data.frame(tips)
dat$new.label <- substr(dat$tips, 1, 20)
dat$new.label <- as.factor(gsub('_$', '', dat$new.label))
tree <- sub.taxa.label(raxml, dat)

ggtree(tree, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=2) +
  geom_rootpoint() # showing node numbers

reroot <- reroot(tree, 196)

# Set-up bootstrap support values
t1 <- ggtree(reroot)
bootstraps_values = t1$data
bootstraps_values = bootstraps_values[!bootstraps_values$isTip,]
bootstraps_values$label = as.numeric(bootstraps_values$label)
bootstraps_values = bootstraps_values[bootstraps_values$label > 0.70,]
bootstraps_values$bs <- with(bootstraps_values, ifelse(bootstraps_values$label>=0.90, 
                                                       "BS>=0.90", ifelse(bootstraps_values$label>=0.70, "BS>=0.70")))
bootstraps_values = na.omit(bootstraps_values)


t2 <- t1 %<+% bootstraps_values + 
  geom_nodepoint(aes(label = bs, color = as.factor(bs)), 
                 size = 1.5, show.legend = FALSE, na.rm = TRUE)
### Associated data sets =======================================

# Metadata
assignment <- t1$data[1:127, 4]
assignment <- assignment[-c(58:60, 125:127),]
assignment$new_label <- substr(assignment$label, 15, 20)
assignment$cluster <- substr(assignment$label, 1, 4)
assignment$taxonomy <- ifelse(grepl("[a-z]", assignment$cluster), "false", 
                              word(assignment$cluster, 1))
assignment$taxonomy <- ifelse(grepl("[A-Z]", assignment$taxonomy), "true",
                              word(assignment$taxonomy, 1))
assignment$cluster <- recode(assignment$cluster,
                             "LSCA" = 1, "Lsca" = 1, "Lgla" = 1.5, "LGLA" = 1.5, "Lmyc" = 2,
                             "LMYC" = 2, "Lhaw" = 2.5, "LHAW" = 2.5, "LAMI" = 3, "LSP1" = 3.5, 
                             "Lspp" = 4, "LINC" = 4, "Lfol" = 4, "Lfra" = 4)

assignment$taxonomy <- as.numeric(as.logical(assignment$taxonomy))
lb = assignment$label
d = data.frame(label=lb, label2 = paste(assignment$new_label), color.cluster = paste(assignment$cluster))

assignment$location <- substr(assignment$label, 6, 9)
assignment$location <- recode(assignment$location, 
                              "CSBL" = 0.2,
                              "CSDT" = 0.2,
                              "CSFL" = 0.2,
                              "CSHO" = 0.2,
                              "CSNW" = 0.2,
                              "GBDA" = 0.1,
                              "GBGD" = 0.1,
                              "GBMY" = 0.1,
                              "GBTI" = 0.1,
                              "GBYO" = 0.1,
                              "HAAC" = 0.3,
                              "HABI" = 0.3,
                              "HAGP" = 0.3,
                              "HAHI" = 0.4,
                              "HAJI" = 0.4,
                              "HALI" = 0.3,
                              "HANH" = 0.3,
                              "HAOA" = 0.3,
                              "HAPA" = 0.4,
                              "ISIU" = 0,
                              "ISNR" = 0,
                              "ISPS" = 0)

assignment$colorlocation <- substr(assignment$label, 1, 1)
assignment$depth <- as.numeric(substr(assignment$label, 11, 13))
assignment$depthcat <- with(assignment, ifelse(depth<31, "Shallow", ifelse(depth<60, "UMCE", ifelse(depth<100, "LMCE", "Deep"))))
assignment$depthcat <- recode(assignment$depthcat, "Shallow" = 8, "UMCE" = 9, "LMCE" = 10, "Deep" = 11)

## DAPC data
DAPC_lepto <- read.csv("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c2 - DAPC/temp-unlinked-snps/agalepto_lepto_denovo_DAPC_posterior_4c2_2.csv",
                header = TRUE, stringsAsFactors = TRUE, check.names = FALSE)
rep1_k6_lepto <- DAPC_lepto[1:120, c(1, 2:7)]
rep1_k7_lepto <- DAPC_lepto[1:120, c(1, 8:14)]
rep1_k8_lepto <- DAPC_lepto[1:120, c(1, 15:22)]
rep1_k9_lepto <- DAPC_lepto[1:120, c(1, 23:31)]
rep1_k10_lepto <- DAPC_lepto[1:120, c(1, 32:41)]
rep1_k11_lepto <- DAPC_lepto[1:120, c(1, 42:52)]
rep1_k12_lepto <- DAPC_lepto[1:120, c(1, 53:64)]
rep1_k13_lepto <- DAPC_lepto[1:120, c(1, 65:77)]
rep1_k14_lepto <- DAPC_lepto[1:120, c(1, 78:91)]

# Format dataframe for plotting
rep1_k6_melt_lepto <- melt(rep1_k6_lepto)
rep1_k7_melt_lepto <- melt(rep1_k7_lepto)
rep1_k8_melt_lepto <- melt(rep1_k8_lepto)
rep1_k9_melt_lepto <- melt(rep1_k9_lepto)
rep1_k10_melt_lepto <- melt(rep1_k10_lepto)
rep1_k11_melt_lepto <- melt(rep1_k11_lepto)
rep1_k12_melt_lepto <- melt(rep1_k12_lepto)
rep1_k13_melt_lepto <- melt(rep1_k13_lepto)
rep1_k14_melt_lepto <- melt(rep1_k14_lepto)

## Admixture data
samples_names<-read.delim("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c3 - Admixture/4c3b - Full dataset SNPS/leptoseris_50MD/agalepto_lepto_4c3.fam",
                          header=FALSE, sep=" ") %>% 
  select(., V2) %>%
  rename(., INDIV=V2)

# read Q file
Qval<-read.table(paste0("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c3 - Admixture/4c3a - Unlinked SNPS/leptoseris-unlinked/agalepto_lepto_4c3_temp.14.Q"))
names(Qval)<-paste0("Cluster", 1:ncol(Qval))

# Format Q file for plotting
# add  sample names  to Qtable
Qval<-cbind(INDIV=samples_names$INDIV, Qval) 

# transform to long format  
Qval_long<- gather(Qval, key=Kgroup, value=Qadmixture, 2:ncol(Qval))

## TSNE data
TSNE <- read.csv("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c1 - PCA-TSNE/agalepto_lepto_tsne_analysis_4c1.csv", header = TRUE)
TSNE_groups <- read.csv("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c2 - DAPC/agalepto_lepto_assignment_4c2.csv", header = FALSE)
TSNE_groups <- rename(TSNE_groups, X = V1)
tsne_merged <- merge(TSNE, TSNE_groups, by = "X")
tsne_merged <- rename(tsne_merged, Individual = X, cluster = V2)

## Densitree data
agalepto_popfile = read.table("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4a - Maximum likelihood inference/agalepto_lepto_popfile_4a.txt",
                              col.names = c("sample", "cluster"), 
                              header = FALSE, stringsAsFactors = FALSE)
clusters <- read.csv("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c2 - DAPC/agalepto_lepto_assignment_4c2.csv", header = FALSE)
clusters <- rename(clusters, sample = V1, cluster2 = V2)
popfile_merged <- agalepto_popfile %>% full_join(clusters, by = "sample", "cluster2") %>%
  select(sample, cluster2, everything())
popfile_merged[is.na(popfile_merged)] <- "OG"

# Color and shape schemes across all plots
color_scheme = c("LG" = "#4C7970", "LG1" = "#5C9992", 
                 "LG2" = "#85ADA8", "LG3" = "#A5CAC6",
                 "LH" = "#6E1210",
                 "LH1" = "#9F2D2D", "LH2" = "#BF7375", 
                 "LH3" = "#B34244", "LH4" = "#370E11", 
                 "LH5" = "#956E71",
                 "LM" = "#FFC037", "LSPP" = "#CCCCCC",
                 "LM1" = "#D9881C", "LM2" = "#FB8500",  
                 "LS" = "#283642", "LS1" = "#4390B9",
                 "LS2" = "#73B5D2", "LS3" = "#a9d6e5",
                 "LS1a" = "#4390B9",  "LS1b" = "#a9d6e5", 
                 "LS1c" = "#64ADCE", "LS0a" = "#121F2B",
                 "LS0b" = "#73B5D2", "LS0c" = "#66A3CC", 
                 "BS>=0.90" = "#231F20", 
                 "BS>=0.70" = "#808285", "OG" = "#CCCCCC",
                 "L" ="#999999", "3" = "#bdbdbd", 
                 "1.5" = "#4C7970", "2.5" = "#6E1210", 
                 "2" = "#FCAF1F", "1" = "#121F2B", 
                 "3.5" = "#bdbdbd", "4" = "#bdbdbd",
                 "8" = "#92B5C8", "9" = "#3E7EA3",
                 "10" = "#24476B", "11" = "#050A0F",
                 "Cluster1" = "#73B5D2", "Cluster2" = "#85ADA8", 
                 "Cluster3" = "#4390B9", "Cluster4" = "#9F2D2D",
                 "Cluster5" = "#BF7375", "Cluster6" = "#5C9992", 
                 "Cluster7" = "#a9d6e5", "Cluster8" = "#B34244", 
                 "Cluster9" = "#64ADCE", "Cluster10" = "#283642", 
                 "Cluster11" = "#4C7970", "Cluster12" = "#D9881C", 
                 "Cluster13" = "#6E1210",  "Cluster14" = "#FFC037")

shape_scheme = c("0" = 13, "1" = 19)
# Individual plots =================================================

tree <- t2 %<+% d + 
  geom_tiplab(aes(label = label2), 
              size = 2, align = TRUE, linesize = 0.3) +
  geom_treescale(x = -0.0001, y = 25, 
                 fontsize = 2.5, offset = 0.6, linesize = 0.5, width = 0.005)
p0 <- facet_plot(tree, panel="Taxonomy", data = assignment, 
                 geom=geom_point2, size = 2, 
                 mapping = aes(x = cluster, color = as.factor(cluster), 
                               shape = as.factor(taxonomy))) +
  scale_color_manual(values = color_scheme, guide = "none")
p1 <- facet_plot(p0, panel = "Depth", data = assignment,
                 geom = geom_point2, size = 2,
                 mapping=aes(x = depthcat, color = as.factor(depthcat))) +
  xlim_expand(c(7.8, 11.2), "Depth")
p2 <- facet_plot(p1, panel="Location", data = assignment, 
                 geom=geom_point2, size = 2, 
                 mapping = aes(x = location, color = as.factor(colorlocation))) +
  xlim_expand(c(0, 0.5), "Location")
p3 <- facet_plot(p2, panel="k=6", data=rep1_k6_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none") 
p4 <- facet_plot(p3, panel="k=7", data=rep1_k7_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none")
p5 <- facet_plot(p4, panel="k=8", data=rep1_k8_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none") 
p6 <- facet_plot(p5, panel="k=9", data=rep1_k9_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none") 
p7 <- facet_plot(p6, panel="k=10", data=rep1_k10_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none")
p8 <- facet_plot(p7, panel="k=11", data=rep1_k11_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none")
p9 <- facet_plot(p8, panel="k=12", data=rep1_k12_melt_lepto, 
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)), 
                 stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none")
p10 <- facet_plot(p9, panel="k=13", data=rep1_k13_melt_lepto, 
                  geom=geom_barh, width = 1.0,
                  mapping = aes(x = value, fill = as.factor(variable)), 
                  stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none")
p11 <- facet_plot(p10, panel="k=14", data=rep1_k14_melt_lepto, 
                  geom=geom_barh, width = 1.0,
                  mapping = aes(x = value, fill = as.factor(variable)), 
                  stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none") 
p12 <- facet_plot(p11, panel="Admixture k=13", data=Qval_long, 
                  geom=geom_barh, width = 1.0,
                  mapping = aes(x = Qadmixture, fill = as.factor(Kgroup)), 
                  stat='identity') +
  scale_fill_manual(values = color_scheme, guide = "none") +
  scale_shape_manual(values = shape_scheme) +
  theme(legend.title = element_blank(), legend.position = "none", strip.background = element_blank())
p13 <- facet_labeller(p12, c(Tree = "Phylogeny")) +
  theme(panel.spacing = unit(0, "lines")) 

# Change the width of the individual facet_plot panels
gt = ggplot_gtable(ggplot_build(p13))
gt$widths[5] = 2*gt$widths[5]
gt$widths[7] = 0.40*gt$widths[7]
gt$widths[9] = 0.25*gt$widths[9]
gt$widths[11] = 0.40*gt$widths[11]
gt$widths[13] = 0.10*gt$widths[13]
gt$widths[15] = 0.10*gt$widths[15]
gt$widths[17] = 0.10*gt$widths[17]
gt$widths[19] = 0.10*gt$widths[19]
gt$widths[21] = 0.10*gt$widths[21]
gt$widths[23] = 0.10*gt$widths[23]
gt$widths[25] = 0.10*gt$widths[25]
gt$widths[27] = 0.10*gt$widths[27]
gt$widths[29] = 0.10*gt$widths[29]
gt$widths[31] = 0.10*gt$widths[31]
grid.draw(gt)
# B: Individual plots =================================================

TSNE_plot <-ggplot(tsne_merged) +
  geom_point(aes(x=X0, y=X1, color=cluster), size = 4) +
  scale_colour_manual(values = color_scheme) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, vjust = 0, hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title.x.bottom = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("TSNE component 1") +
  ylab("TSNE component 2") +
  ggtitle("TSNE k=6")

# C: Individual plots =================================================

plot_density_tree <- function(tree_data, tip_order, popfile, offset){
  ggdensitree(tree_data, alpha = .3,
              layout = 'slanted', 
              tip.order = rev(tip_order)) %<+% popfile +
    geom_tippoint(aes(color = factor(cluster2)),  size = 2) +
    scale_color_manual(values = color_scheme) +
    theme(legend.position='none')
}

tiporder <- function(tree_data){
  with(subset(fortify(tree_data), isTip), label[order(y, decreasing=T)])
}

lepto_tetrad_trees <- read.newick("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4b - SNPs-based species tree/agalepto_lepto_tetrad_4b.tree.boots")
lepto_tetrad_cons_tree <- read.tree("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4b - SNPs-based species tree/agalepto_lepto_tetrad_4b.tree.cons.tree")

ggtree(lepto_tetrad_cons_tree, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=2) +
  geom_rootpoint() # showing node numbers

lepto_tetrad_cons_tree_rooted <- root(lepto_tetrad_cons_tree, node = 157)

lepto_tetrad_cons_tree_flipped <- ggtree(lepto_tetrad_cons_tree_rooted)
lepto_tetrad_cons_tree_flipped_1 <- flip(lepto_tetrad_cons_tree_flipped, 174, 169)                             
lepto_tetrad_cons_tree_flipped_2 <- as.phylo(lepto_tetrad_cons_tree_flipped_1)


ggtree(lepto_tetrad_cons_tree_flipped_1, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=2) +
  geom_rootpoint() # showing node numbers

lepto_tetrad_cons_tree_flipped_1+
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=2) +
  geom_rootpoint() # showing node numbers

plot_tetrad_dens_tree <- plot_density_tree(lepto_tetrad_trees, 
                                           tiporder(lepto_tetrad_cons_tree_flipped_2),
                                           popfile_merged, 0)
# Combine plots
B_C <- plot_grid(TSNE_plot, plot_tetrad_dens_tree, ncol=1, nrow = 2, 
                 labels=c("B", "C"))
A <- plot_grid(gt, B_C, labels = c("A", ""),
               rel_widths = c(1,0.3),
               rel_heights = c(1, 0.1))


# Export plots

ggsave(file="leptoseris_draft_fig2.png", height = 8.27, width = 15, dpi = 300)
plot(A)
dev.off()

svglite(file="leptoseris_draft_fig2.svg", width = 15, height = 8.27, pointsize = 12)
plot(A)
dev.off()
