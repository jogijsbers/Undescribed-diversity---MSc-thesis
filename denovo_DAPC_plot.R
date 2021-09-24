# DAPC from VCF file
# Johanna Gijsbers, Feb 2021

## Dependencies ==============================================================
library("adegenet")
library("vcfR")
library("ggplot2")
library("dplyr")
library("svglite")

## Main code =================================================================
## Import with vcfR
coral.genlight <- vcfR2genlight(read.vcfR("/Volumes/Johanna/Agalepto/Agalepto_manuscript/Agalepto_rad/4 - Leptoseris/4c - Unsupervised clustering analyses/4c2 - DAPC/agalepto_lepto_singlesnp1_4c2.vcf"))
# 4 clusters

grp <- find.clusters(coral.genlight, stat = "BIC", 
                     n.iter = 100000, n.pca = 120, 
                     choose.n.clust = 4,
                     n.start = 1000, set.seed(2))

pop <- grp$grp
max_PCAs <- as.integer(length(pop) / 3) # as <= N/3 advised
dapc.pop <- dapc(coral.genlight, pop, n.pca = max_PCAs, n.da = 6)
optimum_score <- optim.a.score(dapc.pop)

dapc.pop <- dapc(coral.genlight, pop, n.pca = optimum_score$best, n.da = 6)

# 5 clusters
grp5 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 5,
                      n.pca = 50, stat = "BIC", n.iter = 100000, 
                      n.start = 1000, set.seed(2))
pop5 <- grp5$grp
dapc5 <- dapc(coral.genlight, pop5, n.pca = 18,
              n.da = 10)

# 6 clusters
grp6 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 6,
                      n.pca = 50, stat = "BIC", n.iter = 100000, 
                      n.start = 1000, set.seed(2))
pop6 <- grp6$grp
dapc6 <- dapc(coral.genlight, pop6, n.pca = 18,
              n.da = 10)

# 7 clusters
grp7 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 7,
                      n.pca = 50, stat = "BIC", n.iter = 100000, 
                      n.start = 1000, set.seed(2))
pop7 <- grp7$grp
dapc7 <- dapc(coral.genlight, pop7, n.pca = 18,
              n.da = 10)

# 8 clusters
grp8 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 8,
                      n.pca = 50, stat = "BIC", n.iter = 100000, 
                      n.start = 1000, set.seed(2))
pop8 <- grp8$grp
dapc8 <- dapc(coral.genlight, pop8, n.pca = 18,
              n.da = 10)

# 9 clusters
grp9 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 9,
                      n.pca = 50, stat = "BIC", n.iter = 100000, 
                      n.start = 1000, set.seed(2))
pop9 <- grp9$grp
dapc9 <- dapc(coral.genlight, pop9, n.pca = 18,
              n.da = 10)

# 10 clusters
grp10 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 10,
                       n.pca = 50, stat = "BIC", n.iter = 100000, 
                       n.start = 1000, set.seed(2))
pop10 <- grp10$grp
dapc10 <- dapc(coral.genlight, pop10, n.pca = 18,
               n.da = 10)

# 11 clusters
grp11 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 11,
                       n.pca = 50, stat = "BIC", n.iter = 100000, 
                       n.start = 1000, set.seed(2))
pop11 <- grp11$grp
dapc11 <- dapc(coral.genlight, pop11, n.pca = 18,
               n.da = 10)

# 12 clusters
grp12 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 12,
                       n.pca = 50, stat = "BIC", n.iter = 100000, 
                       n.start = 1000, set.seed(2))
pop12 <- grp12$grp
dapc12 <- dapc(coral.genlight, pop12, n.pca = 18,
               n.da = 10)
# 13 clusters
grp13 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 13,
                       n.pca = 50, stat = "BIC", n.iter = 100000, 
                       n.start = 1000, set.seed(2))
pop13 <- grp13$grp
dapc13 <- dapc(coral.genlight, pop13, n.pca = 18,
               n.da = 10)

# 11 clusters
grp11 <- find.clusters(coral.genlight, max.n.clust = 20, choose.n.clust = 11,
                       n.pca = 50, stat = "BIC", n.iter = 100000, 
                       n.start = 1000, set.seed(2))
pop11 <- grp11$grp
dapc11 <- dapc(coral.genlight, pop11, n.pca = 18,
               n.da = 10)


# 15 clusters
grp15 <- find.clusters(coral.genlight, max.n.clust = 20, n.pca = 50, stat = "BIC", n.iter = 100000, n.start = 1000, set.seed(2))
pop15 <- grp15$grp
dapc15 <- dapc(coral.genlight, pop15)

# 16 clusters
grp16 <- find.clusters(coral.genlight, max.n.clust = 20, n.pca = 50, stat = "BIC", n.iter = 100000, n.start = 1000, set.seed(2))
pop16 <- grp16$grp
dapc16 <- dapc(coral.genlight, pop16)
## Plots

# 5 clusters
png("agalepto_aga_denovo_DAPC_5c.png")
scatter(dapc4, label.inds = list(air = 1, pch = NA),
        posi.da="bottomleft", bg="white", pch=16, cstar=1,
        col = c("1" = "darkorange", "2" = "darkslategray", 
                "3" = "mediumseagreen", "4" = "darkred", 
                "5" = "lightslategray"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topleft")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c.png")
compoplot(dapc, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "darkorange", "2" = "darkslategray", 
                  "3" = "mediumseagreen", "4" = "darkred", 
                  "5" = "lightslategray"))
dev.off()

# 6 clusters
png("agalepto_aga_denovo_DAPC_5c_1.png")
scatter(dapc6, label.inds = list(air = 1, pch = NA),
        posi.da="bottomleft", bg="white", pch=16, cstar=1,
        col = c("1" = "lightslategray", "2" = "darkslategray", 
                "3" = "#589897", "4" = "darkorange", 
                "5" = "darkred", "6" = "mediumseagreen"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topleft")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_1.png")
compoplot(dapc6, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightslategray", "2" = "darkslategray", 
                  "3" = "#589897", "4" = "darkorange", 
                  "5" = "darkred", "6" = "mediumseagreen"))
dev.off()


# 7 clusters
png("agalepto_aga_denovo_DAPC_5c_2.png")
scatter(dapc7, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "darkslategray", "2" = "darkorange", 
                "3" = "darkred", "4" = "mediumseagreen", 
                "5" = "#589897", "6" = "lightslategray", "7" = "gold"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_2.png")
compoplot(dapc7, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "darkslategray", "2" = "darkorange", 
                  "3" = "darkred", "4" = "mediumseagreen", 
                  "5" = "#589897", "6" = "lightslategray", "7" = "gold"))
dev.off()

# 8 clusters
png("agalepto_aga_denovo_DAPC_5c_3.png")
scatter(dapc8, label.inds = list(air = 1, pch = NA),
        posi.da="bottomleft", bg="white", pch=16, cstar=1,
        col = c("1" = "lightslategray", "2" = "mediumseagreen", 
                "3" = "darkorange", "4" = "#589897", 
                "5" = "gold", "6" = "darkslategray", 
                "7" = "cadetblue", "8" = "darkred"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topleft")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_3.png")
compoplot(dapc8, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightslategray", "2" = "mediumseagreen", 
                  "3" = "darkorange", "4" = "#589897", 
                  "5" = "gold", "6" = "darkslategray", 
                  "7" = "cadetblue", "8" = "darkred"))
dev.off()

# 9 clusters

png("agalepto_aga_denovo_DAPC_5c_4.png")
scatter(dapc9, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "cadetblue", "2" = "darkred", 
                "3" = "#589897", "4" = "lightseagreen", 
                "5" = "mediumseagreen", "6" = "darkslategray", 
                "7" = "gold", "8" = "darkorange", 
                "9" = "lightslategray"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_4.png")
compoplot(dapc9, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "cadetblue", "2" = "darkred", 
                  "3" = "#589897", "4" = "lightseagreen", 
                  "5" = "mediumseagreen", "6" = "darkslategray", 
                  "7" = "gold", "8" = "darkorange", 
                  "9" = "lightslategray"))
dev.off()

# 10 clusters

png("agalepto_aga_denovo_DAPC_5c_5.png")
scatter(dapc10, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "cadetblue", "2" = "gold", 
                "3" = "mediumseagreen", "4" = "darkseagreen", 
                "5" = "#589897", "6" = "darkslategray", 
                "7" = "lightslategray", "8" = "darkorange", 
                "9" = "lightseagreen", "10" = "darkred"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_5.png")
compoplot(dapc10, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "cadetblue", "2" = "gold", 
                  "3" = "mediumseagreen", "4" = "darkseagreen", 
                  "5" = "#589897", "6" = "darkslategray", 
                  "7" = "lightslategray", "8" = "darkorange", 
                  "9" = "lightseagreen", "10" = "darkred"))
dev.off()

# 11 clusters

png("agalepto_aga_denovo_DAPC_5c_6.png")
scatter(dapc11, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "mediumseagreen", "2" = "darkseagreen", 
                "3" = "darkred", "4" = "lightslategray", 
                "5" = "darkorange", "6" = "cadetblue", 
                "7" = "gold", "8" = "#589897", 
                "9" = "darkslategray", "10" = "midnightblue", 
                "11" = "indianred"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_6.png")
compoplot(dapc11, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "mediumseagreen", "2" = "darkseagreen", 
                  "3" = "darkred", "4" = "lightslategray", 
                  "5" = "darkorange", "6" = "cadetblue", 
                  "7" = "gold", "8" = "#589897", 
                  "9" = "darkslategray", "10" = "midnightblue", 
                  "11" = "indianred"))
dev.off()

# 12 clusters

png("agalepto_aga_denovo_DAPC_5c_7.png")
scatter(dapc12, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "indianred", "2" = "cadetblue", 
                "3" = "#589897", "4" = "darkorange", 
                "5" = "mediumseagreen", "6" = "darkslategray", 
                "7" = "lightseagreen", "8" = "midnightblue", 
                "9" = "gold", "10" = "lightslategray", 
                "11" = "darkred", "12" = "aquamarine"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_7.png")
compoplot(dapc12, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "indianred", "2" = "cadetblue", 
                  "3" = "#589897", "4" = "darkorange", 
                  "5" = "mediumseagreen", "6" = "darkslategray", 
                  "7" = "lightseagreen", "8" = "midnightblue", 
                  "9" = "gold", "10" = "lightslategray", 
                  "11" = "darkred", "12" = "aquamarine"))
dev.off()

# 13 clusters

png("agalepto_aga_denovo_DAPC_5c_8.png")
scatter(dapc13, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "lightslategray", "2" = "darkslategray", 
                "3" = "midnightblue", "4" = "lightseagreen", 
                "5" = "darkseagreen", "6" = "aquamarine", 
                "7" = "cadetblue", "8" = "indianred", 
                "9" = "darkred", "10" = "gold", 
                "11" = "mediumseagreen", "12" = "#589897",
                "13" = "darkorange"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_8.png")
compoplot(dapc13, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightslategray", "2" = "darkslategray", 
                  "3" = "midnightblue", "4" = "lightseagreen", 
                  "5" = "darkseagreen", "6" = "aquamarine", 
                  "7" = "cadetblue", "8" = "indianred", 
                  "9" = "darkred", "10" = "gold", 
                  "11" = "mediumseagreen", "12" = "#589897",
                  "13" = "darkorange"))
dev.off()

# 14 clusters

png("agalepto_aga_denovo_DAPC_5c_9.png")
scatter(dapc14, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "lightseagreen", "2" = "lightslategray", 
                "3" = "#589897", "4" = "#adb5bd", 
                "5" = "cadetblue", "6" = "darkslategray", 
                "7" = "mediumseagreen", "8" = "midnightblue", 
                "9" = "indianred", "10" = "darkorange", 
                "11" = "darkred", "12" = "aquamarine",
                "13" = "darkseagreen", "14" = "gold"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_9.png")
compoplot(dapc14, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightseagreen", "2" = "lightslategray", 
                  "3" = "#589897", "4" = "#adb5bd", 
                  "5" = "cadetblue", "6" = "darkslategray", 
                  "7" = "mediumseagreen", "8" = "midnightblue", 
                  "9" = "indianred", "10" = "darkorange", 
                  "11" = "darkred", "12" = "aquamarine",
                  "13" = "darkseagreen", "14" = "gold"))
dev.off()

# 15 clusters

png("agalepto_aga_denovo_DAPC_5c_10.png")
scatter(dapc15, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "lightseagreen", "2" = "lightslategray", 
                "3" = "#589897", "4" = "#adb5bd", 
                "5" = "cadetblue", "6" = "darkslategray", 
                "7" = "mediumseagreen", "8" = "midnightblue", 
                "9" = "indianred", "10" = "darkorange", 
                "11" = "darkred", "12" = "aquamarine",
                "13" = "darkseagreen", "14" = "gold", "15"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_10.png")
compoplot(dapc15, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightseagreen", "2" = "lightslategray", 
                  "3" = "#589897", "4" = "#adb5bd", 
                  "5" = "cadetblue", "6" = "darkslategray", 
                  "7" = "mediumseagreen", "8" = "midnightblue", 
                  "9" = "indianred", "10" = "darkorange", 
                  "11" = "darkred", "12" = "aquamarine",
                  "13" = "darkseagreen", "14" = "gold", "15" = ""))
dev.off()

# 16 clusters

png("agalepto_aga_denovo_DAPC_5c_11.png")
scatter(dapc16, label.inds = list(air = 1, pch = NA),
        posi.da="bottomright", bg="white", pch=16, cstar=1,
        col = c("1" = "lightseagreen", "2" = "lightslategray", 
                "3" = "#589897", "4" = "#adb5bd", 
                "5" = "cadetblue", "6" = "darkslategray", 
                "7" = "mediumseagreen", "8" = "midnightblue", 
                "9" = "indianred", "10" = "darkorange", 
                "11" = "darkred", "12" = "aquamarine",
                "13" = "darkseagreen", "14" = "gold"),
        scree.pca=TRUE, posi.pca="bottomleft", ratio.pca = 1, inset.pca = 1,
        cell = 2, legend = TRUE, inset.solid = 0.8, solid = 0.9, clab = 0,
        cleg = 1, posi.leg = "topright")
dev.off()

png("agalepto_aga_denovo_DAPC_comp_5c_11.png")
compoplot(dapc16, posi="topright", lab= rownames(coral.genlight),
          ncol=4, cleg = 0.5, show.lab = TRUE,
          col = c("1" = "lightseagreen", "2" = "lightslategray", 
                  "3" = "#589897", "4" = "#adb5bd", 
                  "5" = "cadetblue", "6" = "darkslategray", 
                  "7" = "mediumseagreen", "8" = "midnightblue", 
                  "9" = "indianred", "10" = "darkorange", 
                  "11" = "darkred", "12" = "aquamarine",
                  "13" = "darkseagreen", "14" = "gold"))
dev.off()

## Set output files
filename_base <- "agalepto_lepto_denovo_DAPC4"
output_indcoord_filename <- paste(filename_base, "_indcoord_4c2.csv", sep = "")
output_posterior_filename <- paste(filename_base, "_posterior_4c2.csv", sep = "")
write.table(dapc.pop$ind.coord, file = output_indcoord_filename, sep = ",")
write.table(dapc.pop$posterior, file = output_posterior_filename, sep = ",")
