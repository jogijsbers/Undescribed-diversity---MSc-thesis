# Author: Johanna Gijsbers Alejandre
# Date created: March, 2020

### Create tangle tree of RAD vs COX1 marker

library("dendextend")
library("ggtree")
library("ggplot2")
library("phytools")
library("ape")
library("ggstance")
library("reshape2")
library("gridExtra")
library("cowplot")
library("dplyr")
library("gtable")
library("grid")
library("ggdendro")
library("dendextend")
library("DECIPHER")
library("svglite")
library("tidyverse")
library("aplot")
library("ggnewscale")
library("phylotools")

# Environment

setwd("~/Dropbox/AgaDiversity/Figures/fig2/tangled_tree/")

# Color schemes across all plots

color_scheme = c("BS>=90" = "#1b9e77", "BS>=70" = "#d95f02",
                 "1" = "#bdbdbd", "1.5" = "#bdbdbd", 
                 "2" = "#bdbdbd", "2.5" = "#5384A2", 
                 "3" = "#5A0C0D", "3.5" = "#6d597a", 
                 "4" = "#D13440", "4.5" = "#43385b", 
                 "5" = "#0D1930", "5.5" = "#737373", 
                 "6" = "#bdbdbd")


# Data

cox_tree <- read.nexus("~/Dropbox/AgaDiversity/Figures/fig2/tangled_tree/agalepto_lepto_COX_raxml_nexus.tre")
rad_tree <- read.nexus("~/Dropbox/AgaDiversity/Figures/fig2/tangled_tree/agalepto_lepto_RAD_raxml_nexus.tre")

## COX-tree setup

tips <- cox_tree$tip.label
dat <- data.frame(tips)
dat$new.label <- substr(dat$tips, 1, 18)
dat$new.label <- as.factor(gsub('_$', '', dat$new.label))
cox <- sub.taxa.label(cox_tree, dat)

ggtree(cox, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=1.5) +
  geom_rootpoint() # showing node numbers

cox <- chronos(cox)


## RAD-tree setup

tips <- rad_tree$tip.label
dat <- data.frame(tips)
dat$new.label <- substr(dat$tips, 1, 18)
dat$new.label <- as.factor(gsub('_$', '', dat$new.label))
tree <- sub.taxa.label(rad_tree, dat)

ggtree(tree, layout="rectangular") + 
  geom_text2(aes(label = node)) + 
  geom_tiplab(align=TRUE, size=1.5) +
  geom_rootpoint() # showing node numbers

rad <- chronos(tree)

## Tanglegram


dend <- as.dendrogram(cox)
dend2 <- as.dendrogram(rad)
dl <- dendlist(dend, dend2)

tanglegram <- tanglegram(dl, 
                         warn = FALSE,
                          highlight_distinct_edges  = FALSE,
                          common_subtrees_color_lines = FALSE,
                          highlight_branches_lwd = FALSE,
                         margin_inner=1,
                         lwd = 1
                         )

dl_1 <- dendlist(
  dend2 %>% 
    set("labels_cex", 0.5) %>%
    set("labels_col", value = c("grey"), k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("#000000"), k = 1),
  dend %>% 
    set("labels_cex", 0.5) %>%
    set("labels_col", value = c("grey"), k=1) %>%
    set("branches_lty", 1))

gt <- dendextend::tanglegram(dl_1, sort = TRUE,
                       common_subtrees_color_lines = FALSE,
                       highlight_distinct_edges  = FALSE,
                       highlight_branches_lwd=FALSE,
                       margin_inner=3, lwd=1, fast = TRUE,
                        edge.lwd = 1,
                       intersecting = TRUE) 





# Export plots

png(file="leptoseris_tangle_tree.png", height = 11.69, width = 15, units = 'in', res = 2000)
dendextend::tanglegram(dl_1, sort = TRUE,
                       common_subtrees_color_lines = FALSE,
                       highlight_distinct_edges  = FALSE,
                       highlight_branches_lwd=FALSE,
                       margin_inner=3, lwd=1, fast = TRUE,
                       edge.lwd = 1,
                       intersecting = TRUE) 
dev.off()

svglite(file="leptoseris_tangle_tree.svg", height = 11.69, width = 15, pointsize = 12)
dendextend::tanglegram(dl_1, sort = TRUE,
                       common_subtrees_color_lines = FALSE,
                       highlight_distinct_edges  = FALSE,
                       highlight_branches_lwd=FALSE,
                       margin_inner=3, lwd=1, fast = TRUE,
                       edge.lwd = 1,
                       intersecting = TRUE)
dev.off()
