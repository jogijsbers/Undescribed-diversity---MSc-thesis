## Author: J. Gijsbers Alejandre
## Date: June 2020

## Plot maps of different geographic regions, depending on the latitude/longitud you provide

## Dependencies #####

library(oceanmap)
library(ggplot2)
library(svglite)

### Add map
####### World ####### 

theme_set(theme_bw())
worldmap <- oceanmap:::.get.worldmap(worldmap)
str(worldmap)  
  
lon <- c(-50, 30)
lat <- c(40, -46)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific and Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

####### Caribbean ####### 

lon <- c(-100, -60)
lat <- c(5, 40)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

###### Indo-Pacific ######

lon <- c(30, -135)
lat <- c(40, -46)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific and Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

##### Zoom-in on map ######

lon <- c(30, 45)
lat <- c(10, 30)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Gulf of Aqaba - Red Sea", 
                border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

lon <- c(140, 170)
lat <- c(-40, -5)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Great Barrier Reef", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

lon <- c(175, -140)
lat <- c(-40, -3)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "United States Minor Outlying Islands", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

lon <- c(-165, -150)
lat <- c(15, 25)
figure(width = 9.75, height = 5.28)
plotmap(xlim = lon, ylim = lat, main = "Hawaiian Archipelago", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)

# Export plots

png(file="agalepto_map_fig1.png", height = 2000, width = 2000, units = 'px', pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific and Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.01, cex.ticks=2.3, las = 0)
dev.off()

svglite(file="agalepto_map_fig1.svg", height = 11.69, width = 15, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific and Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

png(file="agalepto_caribbean_fig1.png", height = 11.69, width = 15, units = 'in', res = 2000)
plotmap(xlim = lon, ylim = lat, main = "Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_caribbean_fig1.svg", height = 11.69, width = 15, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Caribbean", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_indo-pacific_supp1_1.svg", width = 9.75, height = 5.28, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

tiff(file="agalepto_indo-pacific_supp1.tiff", width = 3543, height = 2598, units = 'px', pointsize = 70)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()
  
png(file="agalepto_indo-pacific_supp1.png", width = 3543, height = 2598, units = 'px', pointsize = 50)
plotmap(xlim = lon, ylim = lat, main = "Indo-Pacific", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

tiff(file="agalepto_RS_supp1.tiff", height = 3543, width = 2598, units = 'px', pointsize = 70)
plotmap(xlim = lon, ylim = lat, main = "Red Sea", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_RS_supp1.svg", height = 9.75, width = 5.28, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Red Sea", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_GBR_supp1.svg", height = 9.75, width = 5.28, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Great Barrier Reef", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_USMOI_supp1.svg", height = 9.75, width = 5.28, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "United States Minor Outlying Islands", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()

svglite(file="agalepto_HA_supp1.svg", width = 9.75, height = 5.28, pointsize = 12)
plotmap(xlim = lon, ylim = lat, main = "Hawaiian Archipelago", 
        border = "grey", grid = FALSE, bwd = 0.5, cex.lab=0.8, cex.ticks=0.8)
dev.off()
