# Author: Johanna Gijsbers Alejandre
# Date created: June, 2020

### Create composition of Leptoseris figure for Agalepto thesis

## Dependencies ==============================================================
library("tidyverse")
library("ggplot2")
library("ggtree")
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
# Tree set-up
raxml <- read.tree("~/Dropbox/AgaDiversity/agalepto_rad/4 - Sequence-based analyses/4b - RAxML/4b2 - Leptoseris/4b2a - Leptoseris-HC/agalepto_lepto_min10_besttree_bs_4b2a.tre")
tips <- raxml$tip.label
dat <- data.frame(tips)
dat$new.label <- substr(dat$tips, 1, 18)
dat$new.label <- as.factor(gsub('_$', '', dat$new.label))
tree <- sub.taxa.label(raxml, dat)

ggtree(tree, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=1.5) +
  geom_rootpoint() # showing node numbers

reroot <- reroot(tree, 212)

ggtree(reroot) + 
  geom_text2(aes(label = node))
# Set-up bootstrap support values
t1 <- ggtree(reroot)
bootstraps_values = t1$data
bootstraps_values = bootstraps_values[!bootstraps_values$isTip,]
bootstraps_values$label = as.numeric(bootstraps_values$label)
bootstraps_values = bootstraps_values[bootstraps_values$label > 70,]
bootstraps_values$bs <- with(bootstraps_values, ifelse(bootstraps_values$label>=90, "BS>=90", ifelse(bootstraps_values$label>=70, "BS>=70")))
bootstraps_values = na.omit(bootstraps_values)


t2 <- t1 %<+% bootstraps_values + 
  geom_nodepoint(aes(label = bs, color = as.factor(bs)), 
                 size = 1.5, show.legend = TRUE)

# Associated data sets
assignment = read.csv("~/Dropbox/AgaDiversity/Figures/fig3 - Leptoseris main figure/lepto_assignment_1.csv", 
                      header = TRUE, stringsAsFactors = TRUE)
assignment$Cluster <- recode(assignment$Cluster,
                             "LS" = 1, "Ls" = 1, "Lg" = 2, "LG" = 2, "Lm" = 3, "LM" = 3, 
                             "Lh" = 4, "LH" = 4, "LA" = 5, "LZ" = 6, "LX" = 7)

assignment$Taxonomy <- as.numeric(assignment$Taxonomy)
lb = assignment$Individual
d = data.frame(label=lb, label2 = paste(assignment$New.label))


assignment$Location <- substr(assignment$Individual, 4, 5)
assignment$Location <- recode(assignment$Location, "IS" = 0, "GB" = 1, "CS" = 2, "HA" = 3)
assignment$colorlocation <- substr(assignment$Individual, 1, 1)
assignment$colorlocation <- recode(assignment$colorlocation, "C" = "A", "L" = "A")
assignment$colorlocation <- as.factor(assignment$colorlocation)

assignment$pop <- as.factor(substr(assignment$Individual, 1, 2))
assignment$depth <- substr(assignment$Individual, 9,11)
assignment$depth <- as.factor(gsub('_$', '', assignment$depth))
assignment$depth <- as.numeric(as.character(assignment$depth))
assignment$depthcat <- with(assignment, ifelse(depth<31, "Shallow", ifelse(depth<60, "UMCE", ifelse(depth<100, "LMCE", "Deep"))))
assignment$depthcat <- recode(assignment$depthcat, "Shallow" = 8, "UMCE" = 9, "LMCE" = 10, "Deep" = 11)

heatmapDepth = data.frame(assignment$Individual, row.names = 1)
heatmapDepth$Depth <- as.numeric(substr(assignment$Individual, 9,11))
heatmapDepth$Depth <- as.factor(gsub('_$', '', assignment$depth))
heatmapDepth$Depth <- as.numeric(as.character(assignment$depth))

heatmapData=read.csv("~/Dropbox/AgaDiversity/Figures/fig3 - Leptoseris main figure/leptoseris_spp_del_models_ranked_2.csv", row.names = 1)
rn <- rownames(heatmapData)
heatmapData <- heatmapData[,c(23:27)]
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn

DAPC = read.csv("~/Dropbox/AgaDiversity/agalepto_rad/3 - SNP-based analyses/3c - DAPC/3c2 - Leptoseris/3c2a - Leptoseris spp./agalepto_lepto_DAPC_posterior_3c2a.csv",
                header = TRUE, stringsAsFactors = TRUE)

dapc_clusters <- melt(DAPC)

# Color schemes across all plots
color_scheme = c("LA" = "#ffffff", "LG1" = "#4C7970", 
                 "LG9" = "#6FBEA8", "LG6" = "#C6E7E2",
                 "LG2" = "#122B27", "LG3" = "#A9DBD3", 
                 "LG8" = "#307368", "LG7" = "#50AE8D",
                 "LG5" = "#7FC8BC","LG4" = "#265650", 
                 "LG" = "#4C7970", "LH" = "#A9504E",
                 "LH1" = "#A9504E", "LH1b" = "#A62B27", 
                 "LH2" = "#BE7774", "LH3" = "#6E1210", 
                 "LH5" = "#D8ACAB", "LM" = "#FCAF1F", 
                 "LM1" = "#FCAF1F", "LM2" = "#D67C28", 
                 "LM3" = "#F3C417", "LO" = "#bdbdbd", 
                 "LS" = "#121F2B", "LS4" = "#2C7CA3", 
                 "LS6" = "#9ECDE6", "LS2c" = "#1F486B", 
                 "LS2b" = "#90A3B6", "LS2" = "#5591B7", 
                 "LS1" = "#121F2B", "LS5" = "#415F79", 
                 "LS3" = "#AFC5D6", "LX" = "#ffffff", 
                 "BS>=90" = "#231F20", "BS>=70" = "#808285", 
                 "A" ="#999999", "5" = "#bdbdbd", 
                 "2" = "#4C7970", "4" = "#6E1210", 
                 "3" = "#FCAF1F", "1" = "#121F2B", 
                 "6" = "#bdbdbd", "7" = "#bdbdbd",
                 "8" = "#92B5C8", "9" = "#3E7EA3",
                 "10" = "#24476B", "11" = "#050A0F")


shape_scheme = c("0" = 13, "1" = 19)
# Individual plots =================================================

tree <- t2 %<+% d + 
  geom_tiplab(aes(label = label2), 
              size = 2, align = TRUE, linesize = 0.3) +
  geom_tippoint(aes(color = ), size=3, alpha=.75) +
  geom_treescale(x = -0.0001, y = 25, 
                 fontsize = 2.5, offset = 0.6, linesize = 0.5, width = 0.005) +
  geom_hilight(node=146, fill="#121F2B", alpha=.8) +
  geom_hilight(node=238, fill="#121F2B", alpha=.8) +
  geom_hilight(node=177, fill="#4C7970", alpha=.8) +
  geom_hilight(node=194, fill="#4C7970", alpha=.8) +
  geom_hilight(node=197, fill="#4C7970", alpha=.8) +
  geom_hilight(node=202, fill="#FCAF1F", alpha=.8) +
  geom_hilight(node=208, fill="#A84F4D", alpha=.8) +
  geom_hilight(node=222, fill="#A84F4D", alpha=.8) +
  geom_hilight(node=228, fill="#A84F4D", alpha=.8) +
  geom_hilight(node=251, fill="#FCAF1F", alpha=.8) +
  scale_colour_discrete(na.translate = F)

p1 <- facet_plot(tree, panel = "Depth", data = assignment,
                 geom = geom_point, size = 2,
                 aes(x = depthcat, color = as.factor(depthcat))) +
  theme(panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = color_scheme, guide = "none") +
  theme(legend.position = c(0.05, 0.8),
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")) +
  labs(fill = "Depth (m)", color = "Support")

p1
  
p2 <- p1 + new_scale_fill()

p3 <- gheatmap(p2, heatmapData, color = "white", offset = 0.01, 
               width = 0.3, font.size = 4, colnames_position = 'top', colnames_offset_y = 1) +
  #scale_y_continuous(expand = c(0.2,-1)) +
  scale_color_manual(values = color_scheme, guide = "none") +
  theme(panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = color_scheme, guide = "none") +
  theme(legend.position = c(0.05, 0.8),
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")) +
  labs(fill = "Depth (m)")
p3

p4 <- facet_plot(p3, panel="Location", data = assignment, 
                 geom=geom_point, size = 2,
                 aes(x = Location, color = as.factor(colorlocation))) +
  theme(panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = color_scheme, guide = "none") +
  theme(legend.position = c(0.05, 0.8),
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")) +
  labs(fill = "Depth (m)")

p5 <- facet_plot(p4, panel="DAPC", data=dapc_clusters, 
                       geom=geom_barh, width = 1.0,
                       mapping = aes(x = value, fill = as.factor(variable)), 
                       stat='identity') +
  theme(panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = color_scheme, guide = "none") +
  theme(legend.position = c(0.05, 0.8),
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")) +
  labs(fill = "Depth (m)")
  
p6 <- facet_plot(p5, panel="Taxonomy", data = assignment, 
                 geom=geom_point2, size = 2, 
                 aes(x = Cluster, color = as.factor(Cluster), shape = as.factor(Taxonomy))) +
  theme(panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = color_scheme, guide = "none") +
  scale_shape_manual(values = shape_scheme) +
  theme(legend.position = c(0.05, 0.8),
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")) +
  labs(fill = "Depth (m)")

p6
p7 <- facet_labeller(p5, c(Tree = "Phylogeny")) + theme(panel.spacing = unit(0.5, "lines")) 
# Change the width of the individual facet_plot panels
gt = ggplot_gtable(ggplot_build(p6))
gt$widths[5] = 5*gt$widths[5]
gt$widths[7] = 0.35*gt$widths[7]
gt$widths[9] = 0.35*gt$widths[9]
gt$widths[11] = 0.50*gt$widths[11]
grid.draw(gt)

# Export plots

png(file="leptoseris_BFD_fig3.png", height = 11.69, width = 15, units = 'in', res = 2000)
plot(gt)
dev.off()

svglite("leptoserisBFDfig3.svg", height = 11.69, width = 15, pointsize = 12)
plot(gt)
dev.off()
 
