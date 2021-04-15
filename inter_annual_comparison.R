### This script is for analyzing percent annual snow cover by landscape feature and plotting these data for all water years


library(sp)
library(sf)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(caret)
library(data.table)
library(stats)
library(fields)
library(hydroTSM)
library(SpatialTools)
library(foreign)
library(rasterVis)
library(RColorBrewer)
library(viridis)


ER<- readOGR("./EastRiver_Project.shp") #shapefile of the study watershed for reference



######### Compare number percent of snow-covered days to topographic features by water year

## Load in the rasters with elevation, slope, and aspect data (slope and aspect data calculated in dem_anaysis.R script)
east_DEM<- raster('./East_DEM.tif')

east_slope<- raster('./east_slope.tif')

east_aspect<- raster('./east_aspect.tif')


####### Elevation analyses 

# Find the min and max elevations

elev_stats<- stats(east_DEM)
elev_min<- elev_stats[[4]]
elev_max<- elev_stats[[8]]

elev_range<- elev_max-elev_min


## Calculate the empirical cumulative distribution of the elevation data

DEM_vec<- as.vector(east_DEM)

x<- sort(DEM_vec)

e_cdf<- 1:length(x)/ length(x)
#plot(x, e_cdf, type = 's')

low_elev_limit<- x[which(e_cdf >= 0.33)[1]]
mid_elev_limit<- x[which(e_cdf >= 0.66)[1]]


# Create the elevation bands with the ecdf values

elev_low<- (east_DEM <= low_elev_limit)
elev_mid<- (east_DEM > low_elev_limit & east_DEM <= mid_elev_limit)
elev_high<- (east_DEM > mid_elev_limit)

elev_bands <- east_DEM
elev_bands[elev_low] <- 1
elev_bands[elev_mid] <- 2
elev_bands[elev_high] <- 3

plot(elev_bands, main = "Elevation Bands")
plot(ER_proj, border = 'black', add = TRUE)

elev_bands_fac<- as.factor(elev_bands)
rat<- levels(elev_bands_fac)[[1]] 
rat[["bands"]] <- c("low", "medium", "high")
levels(elev_bands_fac)<- rat

levelplot(elev_bands_fac, col.regions = rev(topo.colors(3)), main = "Elevation Bands") + layer(sp.polygons(ER_proj))

#### Percent snow days by elevation band for all water years

### Low elevation

low_elev_perc_2008 <- corrected_2008
masked <- Which(elev_low == 0, cells=TRUE) ## identify all pixels that are not located in the low elevation band
low_elev_perc_2008[masked] <- NA ## mask out these pixels so only low elevation pixels have data

plot(low_elev_perc_2008, main = "2008 Perc snow at Low Elevation")
plot(ER_proj, border = 'black', add = TRUE)

low_elev_perc_2010 <- corrected_2010
masked <- Which(elev_low == 0, cells=TRUE)
low_elev_perc_2010[masked] <- NA

plot(low_elev_perc_2010, main = "2010 Perc snow at Low Elevation")
plot(ER_proj, border = 'black', add = TRUE)

low_elev_perc_2012 <- corrected_2012
masked <- Which(elev_low == 0, cells=TRUE)
low_elev_perc_2012[masked] <- NA

plot(low_elev_perc_2012, main = "2012 Perc snow at Low Elevation")
plot(ER_proj, border = 'black', add = TRUE)

### Medium elevation

mid_elev_perc_2008 <- corrected_2008
masked <- Which(elev_mid == 0, cells=TRUE) ## identify all pixels that are not located in the mid elevation band
mid_elev_perc_2008[masked] <- NA ## mask out these pixels so only mid elevation pixels have data

plot(mid_elev_perc_2008, main = "2008 Perc snow at Mid Elevation")
plot(ER_proj, border = 'black', add = TRUE)

mid_elev_perc_2010 <- corrected_2010
masked <- Which(elev_mid == 0, cells=TRUE)
mid_elev_perc_2010[masked] <- NA

plot(mid_elev_perc_2010, main = "2010 Perc snow at Mid Elevation")
plot(ER_proj, border = 'black', add = TRUE)

mid_elev_perc_2012 <- corrected_2012
masked <- Which(elev_mid == 0, cells=TRUE)
mid_elev_perc_2012[masked] <- NA

plot(mid_elev_perc_2012, main = "2012 Perc snow at Mid Elevation")
plot(ER_proj, border = 'black', add = TRUE)


### High elevation

high_elev_perc_2008 <- corrected_2008
masked <- Which(elev_high == 0, cells=TRUE) ## identify all pixels that are not located in the high elevation band
high_elev_perc_2008[masked] <- NA ## mask out these pixels so only high elevation pixels have data 

plot(high_elev_perc_2008, main = "2008 Perc snow at High Elevation")
plot(ER_proj, border = 'black', add = TRUE)

high_elev_perc_2010 <- corrected_2010
masked <- Which(elev_high == 0, cells=TRUE)
high_elev_perc_2010[masked] <- NA

plot(high_elev_perc_2010, main = "2010 Perc snow at High Elevation")
plot(ER_proj, border = 'black', add = TRUE)

high_elev_perc_2012 <- corrected_2012
masked <- Which(elev_high == 0, cells=TRUE)
high_elev_perc_2012[masked] <- NA

plot(high_elev_perc_2012, main = "2012 Perc snow at High Elevation")
plot(ER_proj, border = 'black', add = TRUE)


### Create graphs to compare average snow percentages by elevation

## Turn all rasters into vectors

low_elev_perc_2008_vec<- as.vector(low_elev_perc_2008)
mid_elev_perc_2008_vec<- as.vector(mid_elev_perc_2008)
high_elev_perc_2008_vec<- as.vector(high_elev_perc_2008)

low_elev_perc_2010_vec<- as.vector(low_elev_perc_2010)
mid_elev_perc_2010_vec<- as.vector(mid_elev_perc_2010)
high_elev_perc_2010_vec<- as.vector(high_elev_perc_2010)

low_elev_perc_2012_vec<- as.vector(low_elev_perc_2012)
mid_elev_perc_2012_vec<- as.vector(mid_elev_perc_2012)
high_elev_perc_2012_vec<- as.vector(high_elev_perc_2012)


## Calculate means and standard deviations

low_elev_2008_mean<- mean(low_elev_perc_2008_vec, na.rm = TRUE)
low_elev_2008_sd<- sd(low_elev_perc_2008_vec, na.rm = TRUE)
mid_elev_2008_mean<- mean(mid_elev_perc_2008_vec, na.rm = TRUE)
mid_elev_2008_sd<- sd(mid_elev_perc_2008_vec, na.rm = TRUE)
high_elev_2008_mean<- mean(high_elev_perc_2008_vec, na.rm = TRUE)
high_elev_2008_sd<- sd(high_elev_perc_2008_vec, na.rm = TRUE)

low_elev_2010_mean<- mean(low_elev_perc_2010_vec, na.rm = TRUE)
low_elev_2010_sd<- sd(low_elev_perc_2010_vec, na.rm = TRUE)
mid_elev_2010_mean<- mean(mid_elev_perc_2010_vec, na.rm = TRUE)
mid_elev_2010_sd<- sd(mid_elev_perc_2010_vec, na.rm = TRUE)
high_elev_2010_mean<- mean(high_elev_perc_2010_vec, na.rm = TRUE)
high_elev_2010_sd<- sd(high_elev_perc_2010_vec, na.rm = TRUE)

low_elev_2012_mean<- mean(low_elev_perc_2012_vec, na.rm = TRUE)
low_elev_2012_sd<- sd(low_elev_perc_2012_vec, na.rm = TRUE)
mid_elev_2012_mean<- mean(mid_elev_perc_2012_vec, na.rm = TRUE)
mid_elev_2012_sd<- sd(mid_elev_perc_2012_vec, na.rm = TRUE)
high_elev_2012_mean<- mean(high_elev_perc_2012_vec, na.rm = TRUE)
high_elev_2012_sd<- sd(high_elev_perc_2012_vec, na.rm = TRUE)

## Create and display the bar plots

## Elevation bar plots with means and standard deviations

labels<- c("2008", '2010', '2012')
means<- c(low_elev_2008_mean, low_elev_2010_mean, low_elev_2012_mean)
sds<- c(low_elev_2008_sd, low_elev_2010_sd, low_elev_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - Low Elev", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
#text(mids, means + 20, paste(round(means, digits = 2)), cex = 1.2)
text(mids, means + 20, paste(round(means)), cex = 1.2)


labels<- c("2008", '2010', '2012')
means<- c(mid_elev_2008_mean, mid_elev_2010_mean, mid_elev_2012_mean)
sds<- c(mid_elev_2008_sd, mid_elev_2010_sd, mid_elev_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - Mid Elev", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 22, paste(round(means)), cex = 1.2)


labels<- c("2008", '2010', '2012')
means<- c(high_elev_2008_mean, high_elev_2010_mean, high_elev_2012_mean)
sds<- c(high_elev_2008_sd, high_elev_2010_sd, high_elev_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - High Elev", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 23, paste(round(means)), cex = 1.2)


################# Slope analyses

slope_stats<- stats(east_slope)
slope_max<- slope_stats[[8]]
slope_min<- slope_stats[[4]]


## Calculate the empirical cumulative distribution of the elevation data

slope_vec<- as.vector(east_slope)

x<- sort(slope_vec)

e_cdf<- 1:length(x)/ length(x)
#plot(x, e_cdf, type = 's')

low_slope_limit<- x[which(e_cdf >= 0.33)[1]]
mid_slope_limit<- x[which(e_cdf >= 0.66)[1]]


# Create the slope bands
slope_low<- (east_slope <= low_slope_limit)
slope_mid<- (east_slope > low_slope_limit & east_slope <= mid_slope_limit)
slope_high<- (east_slope > mid_slope_limit)

slope_bands <- east_slope
slope_bands[slope_low] <- 1
slope_bands[slope_mid] <- 2
slope_bands[slope_high] <- 3

plot(slope_bands, main = "Slope Bands")
plot(ER_proj, border = 'black', add = TRUE)

slope_bands_fac<- as.factor(slope_bands)
rat<- levels(slope_bands_fac)[[1]] 
rat[["bands"]] <- c("low", "medium", "high")
levels(slope_bands_fac)<- rat

levelplot(slope_bands_fac, col.regions = (terrain.colors(3)), main = "Slope Bands") + layer(sp.polygons(ER_proj))


#### Percent snow days by slope band for all water years

# For low slopes

low_slope_perc_2008 <- corrected_2008
masked <- Which(slope_low == 0, cells=TRUE) ## identify all pixels that are not located in the low slope band
low_slope_perc_2008[masked] <- NA ## mask out these pixels so only low slope pixels have data 

## Need to crop the raster due to edge effects from slope calcualtion
low_slope_perc_2008_crop<- crop(low_slope_perc_2008, extent(low_slope_perc_2008, 2, 1590, 2, 1805))

plot(low_slope_perc_2008_crop, main = "2008 Perc snow at Low Slopes")
plot(ER_proj, border = 'black', add = TRUE)


low_slope_perc_2010 <- corrected_2010
masked <- Which(slope_low == 0, cells=TRUE)
low_slope_perc_2010[masked] <- NA

low_slope_perc_2010_crop<- crop(low_slope_perc_2010, extent(low_slope_perc_2010, 2, 1590, 2, 1805))

plot(low_slope_perc_2010_crop, main = "2010 Perc snow at Low Slopes")
plot(ER_proj, border = 'black', add = TRUE)


low_slope_perc_2012 <- corrected_2012
masked <- Which(slope_low == 0, cells=TRUE)
low_slope_perc_2012[masked] <- NA

low_slope_perc_2012_crop<- crop(low_slope_perc_2012, extent(low_slope_perc_2012, 2, 1590, 2, 1805))

plot(low_slope_perc_2012_crop, main = "2012 Perc snow at Low Slopes")
plot(ER_proj, border = 'black', add = TRUE)

# For med slopes

mid_slope_perc_2008 <- corrected_2008
masked <- Which(slope_mid == 0, cells=TRUE) ## identify all pixels that are not located in the mid slope band
mid_slope_perc_2008[masked] <- NA ## mask out these pixels so only mid slope pixels have data

## Need to crop the raster due to edge effects from slope calcualtion
mid_slope_perc_2008_crop<- crop(mid_slope_perc_2008, extent(mid_slope_perc_2008, 2, 1590, 2, 1805))

plot(mid_slope_perc_2008_crop, main = "2008 Perc snow at Mid Slopes")
plot(ER_proj, border = 'black', add = TRUE)


mid_slope_perc_2010 <- corrected_2010
masked <- Which(slope_mid == 0, cells=TRUE)
mid_slope_perc_2010[masked] <- NA

mid_slope_perc_2010_crop<- crop(mid_slope_perc_2010, extent(mid_slope_perc_2010, 2, 1590, 2, 1805))

plot(mid_slope_perc_2010_crop, main = "2010 Perc snow at Mid Slopes")
plot(ER_proj, border = 'black', add = TRUE)


mid_slope_perc_2012 <- corrected_2012
masked <- Which(slope_mid == 0, cells=TRUE)
mid_slope_perc_2012[masked] <- NA

mid_slope_perc_2012_crop<- crop(mid_slope_perc_2012, extent(mid_slope_perc_2012, 2, 1590, 2, 1805))

plot(mid_slope_perc_2012_crop, main = "2012 Perc snow at Mid Slopes")
plot(ER_proj, border = 'black', add = TRUE)

 # For high slopes

high_slope_perc_2008 <- corrected_2008
masked <- Which(slope_high == 0, cells=TRUE) ## identify all pixels that are not located in the high slope band
high_slope_perc_2008[masked] <- NA ## mask out these pixels so only high slope pixels have data

## Need to crop the raster due to edge effects from slope calcualtion
high_slope_perc_2008_crop<- crop(high_slope_perc_2008, extent(high_slope_perc_2008, 2, 1590, 2, 1805))

plot(high_slope_perc_2008_crop, main = "2008 Perc snow at High Slopes")
plot(ER_proj, border = 'black', add = TRUE)


high_slope_perc_2010 <- corrected_2010
masked <- Which(slope_high == 0, cells=TRUE)
high_slope_perc_2010[masked] <- NA

high_slope_perc_2010_crop<- crop(high_slope_perc_2010, extent(high_slope_perc_2010, 2, 1590, 2, 1805))

plot(high_slope_perc_2010_crop, main = "2010 Perc snow at High Slopes")
plot(ER_proj, border = 'black', add = TRUE)


high_slope_perc_2012 <- corrected_2012
masked <- Which(slope_high == 0, cells=TRUE)
high_slope_perc_2012[masked] <- NA

high_slope_perc_2012_crop<- crop(high_slope_perc_2012, extent(high_slope_perc_2012, 2, 1590, 2, 1805))

plot(high_slope_perc_2012_crop, main = "2012 Perc snow at High Slopes")
plot(ER_proj, border = 'black', add = TRUE)


### Create graphs to compare average snow percentages by slope

## Turn all rasters into vectors

low_slope_perc_2008_vec<- as.vector(low_slope_perc_2008_crop)
low_slope_perc_2010_vec<- as.vector(low_slope_perc_2010_crop)
low_slope_perc_2012_vec<- as.vector(low_slope_perc_2012_crop)

mid_slope_perc_2008_vec<- as.vector(mid_slope_perc_2008_crop)
mid_slope_perc_2010_vec<- as.vector(mid_slope_perc_2010_crop)
mid_slope_perc_2012_vec<- as.vector(mid_slope_perc_2012_crop)

high_slope_perc_2008_vec<- as.vector(high_slope_perc_2008_crop)
high_slope_perc_2010_vec<- as.vector(high_slope_perc_2010_crop)
high_slope_perc_2012_vec<- as.vector(high_slope_perc_2012_crop)


## Calculate means and standard deviations

low_slope_2008_mean<- mean(low_slope_perc_2008_vec, na.rm = TRUE)
low_slope_2008_sd<- sd(low_slope_perc_2008_vec, na.rm = TRUE)
low_slope_2010_mean<- mean(low_slope_perc_2010_vec, na.rm = TRUE)
low_slope_2010_sd<- sd(low_slope_perc_2010_vec, na.rm = TRUE)
low_slope_2012_mean<- mean(low_slope_perc_2012_vec, na.rm = TRUE)
low_slope_2012_sd<- sd(low_slope_perc_2012_vec, na.rm = TRUE)

mid_slope_2008_mean<- mean(mid_slope_perc_2008_vec, na.rm = TRUE)
mid_slope_2008_sd<- sd(mid_slope_perc_2008_vec, na.rm = TRUE)
mid_slope_2010_mean<- mean(mid_slope_perc_2010_vec, na.rm = TRUE)
mid_slope_2010_sd<- sd(mid_slope_perc_2010_vec, na.rm = TRUE)
mid_slope_2012_mean<- mean(mid_slope_perc_2012_vec, na.rm = TRUE)
mid_slope_2012_sd<- sd(mid_slope_perc_2012_vec, na.rm = TRUE)

high_slope_2008_mean<- mean(high_slope_perc_2008_vec, na.rm = TRUE)
high_slope_2008_sd<- sd(high_slope_perc_2008_vec, na.rm = TRUE)
high_slope_2010_mean<- mean(high_slope_perc_2010_vec, na.rm = TRUE)
high_slope_2010_sd<- sd(high_slope_perc_2010_vec, na.rm = TRUE)
high_slope_2012_mean<- mean(high_slope_perc_2012_vec, na.rm = TRUE)
high_slope_2012_sd<- sd(high_slope_perc_2012_vec, na.rm = TRUE)


### Create a bar graph with slope means and standard deviations

labels<- c("2008", '2010', '2012')
means<- c(low_slope_2008_mean, low_slope_2010_mean, low_slope_2012_mean)
sds<- c(low_slope_2008_sd, low_slope_2010_sd, low_slope_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - Low Slopes", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 23, paste(round(means)), cex = 1.2)


labels<- c("2008", '2010', '2012')
means<- c(mid_slope_2008_mean, mid_slope_2010_mean, mid_slope_2012_mean)
sds<- c(mid_slope_2008_sd, mid_slope_2010_sd, mid_slope_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - Mid Slopes", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


labels<- c("2008", '2010', '2012')
means<- c(high_slope_2008_mean, high_slope_2010_mean, high_slope_2012_mean)
sds<- c(high_slope_2008_sd, high_slope_2010_sd, high_slope_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow Cover - High Slopes", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim=range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


########### Analyses by aspect    


## Assign directions to variables

aspect_stats<- stats(east_aspect)


north<- (east_aspect >=0 & east_aspect <= 45 | east_aspect > 315)
east<- (east_aspect > 45 & east_aspect <= 135)
south<- (east_aspect > 135 & east_aspect <= 225)
west<- (east_aspect > 225 & east_aspect <= 315)


aspect_cat <- east_aspect
aspect_cat[north] <- 1
aspect_cat[east] <- 2
aspect_cat[south] <- 3
aspect_cat[west] <- 4


aspect_cat<- ratify(aspect_cat)

rat_aspect<- levels(aspect_cat)[[1]]
rat_aspect$direction<- c('North', 'East', 'South', 'West')
levels(aspect_cat)<- rat_aspect

levelplot(aspect_cat, att = 'direction', col.regions = terrain.colors(4), main = 'Aspect') + layer(sp.polygons(ER_proj))


# For north facing slopes

north_aspect_perc_2008 <- corrected_2008
masked <- Which(north == 0, cells=TRUE) ## identify all pixels that are not located on north aspects
north_aspect_perc_2008[masked] <- NA ## mask out these pixels so only north aspect pixels have data

## Need to crop the raster due to edge effects from aspect calcualtion
north_aspect_perc_2008_crop<- crop(north_aspect_perc_2008, extent(north_aspect_perc_2008, 2, 1590, 2, 1805))

plot(north_aspect_perc_2008_crop, main = "2008 Perc snow at north facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


north_aspect_perc_2010 <- corrected_2010
masked <- Which(north == 0, cells=TRUE)
north_aspect_perc_2010[masked] <- NA

north_aspect_perc_2010_crop<- crop(north_aspect_perc_2010, extent(north_aspect_perc_2010, 2, 1590, 2, 1805))

plot(north_aspect_perc_2010_crop, main = "2010 Perc snow at north facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


north_aspect_perc_2012 <- corrected_2012
masked <- Which(north == 0, cells=TRUE)
north_aspect_perc_2012[masked] <- NA

north_aspect_perc_2012_crop<- crop(north_aspect_perc_2012, extent(north_aspect_perc_2012, 2, 1590, 2, 1805))

plot(north_aspect_perc_2012_crop, main = "2012 Perc snow at north facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


# For east facing slopes

east_aspect_perc_2008 <- corrected_2008 
masked <- Which(east == 0, cells=TRUE) ## identify all pixels that are not located on east aspects
east_aspect_perc_2008[masked] <- NA ## mask out these pixels so only east aspect pixels have data

## Need to crop the raster due to edge effects from aspect calcualtion
east_aspect_perc_2008_crop<- crop(east_aspect_perc_2008, extent(east_aspect_perc_2008, 2, 1590, 2, 1805))

plot(east_aspect_perc_2008_crop, main = "2008 Perc snow at east facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


east_aspect_perc_2010 <- corrected_2010
masked <- Which(east == 0, cells=TRUE)
east_aspect_perc_2010[masked] <- NA

east_aspect_perc_2010_crop<- crop(east_aspect_perc_2010, extent(east_aspect_perc_2010, 2, 1590, 2, 1805))

plot(east_aspect_perc_2010_crop, main = "2010 Perc snow at east facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


east_aspect_perc_2012 <- corrected_2012
masked <- Which(east == 0, cells=TRUE)
east_aspect_perc_2012[masked] <- NA

east_aspect_perc_2012_crop<- crop(east_aspect_perc_2012, extent(east_aspect_perc_2012, 2, 1590, 2, 1805))

plot(east_aspect_perc_2012_crop, main = "2012 Perc snow at east facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


# For south facing slopes

south_aspect_perc_2008 <- corrected_2008
masked <- Which(south == 0, cells=TRUE) ## identify all pixels that are not located on south aspects
south_aspect_perc_2008[masked] <- NA ## mask out these pixels so only south aspect pixels have data

## Need to crop the raster due to edge effects from aspect calcualtion
south_aspect_perc_2008_crop<- crop(south_aspect_perc_2008, extent(south_aspect_perc_2008, 2, 1590, 2, 1805))

plot(south_aspect_perc_2008_crop, main = "2008 Perc snow at south facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


south_aspect_perc_2010 <- corrected_2010
masked <- Which(south == 0, cells=TRUE)
south_aspect_perc_2010[masked] <- NA

south_aspect_perc_2010_crop<- crop(south_aspect_perc_2010, extent(south_aspect_perc_2010, 2, 1590, 2, 1805))

plot(south_aspect_perc_2010_crop, main = "2010 Perc snow at south facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


south_aspect_perc_2012 <- corrected_2012
masked <- Which(south == 0, cells=TRUE)
south_aspect_perc_2012[masked] <- NA

south_aspect_perc_2012_crop<- crop(south_aspect_perc_2012, extent(south_aspect_perc_2012, 2, 1590, 2, 1805))

plot(south_aspect_perc_2012_crop, main = "2012 Perc snow at south facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


# For west facing slopes

west_aspect_perc_2008 <- corrected_2008
masked <- Which(west == 0, cells=TRUE) ## identify all pixels that are not located on west aspects
west_aspect_perc_2008[masked] <- NA ## mask out these pixels so only west aspect pixels have data

## Need to crop the raster due to edge effects from aspect calcualtion
west_aspect_perc_2008_crop<- crop(west_aspect_perc_2008, extent(west_aspect_perc_2008, 2, 1590, 2, 1805))

plot(west_aspect_perc_2008_crop, main = "2008 Perc snow at west facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


west_aspect_perc_2010 <- corrected_2010
masked <- Which(west == 0, cells=TRUE)
west_aspect_perc_2010[masked] <- NA

west_aspect_perc_2010_crop<- crop(west_aspect_perc_2010, extent(west_aspect_perc_2010, 2, 1590, 2, 1805))

plot(west_aspect_perc_2010_crop, main = "2010 Perc snow at west facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


west_aspect_perc_2012 <- corrected_2012
masked <- Which(west == 0, cells=TRUE)
west_aspect_perc_2012[masked] <- NA

west_aspect_perc_2012_crop<- crop(west_aspect_perc_2012, extent(west_aspect_perc_2012, 2, 1590, 2, 1805))

plot(west_aspect_perc_2012_crop, main = "2012 Perc snow at west facing slopes")
plot(ER_proj, border = 'black', add = TRUE)


## Turn all rasters into vectors

north_aspect_perc_2008_vec<- as.vector(north_aspect_perc_2008_crop)
north_aspect_perc_2010_vec<- as.vector(north_aspect_perc_2010_crop)
north_aspect_perc_2012_vec<- as.vector(north_aspect_perc_2012_crop)

east_aspect_perc_2008_vec<- as.vector(east_aspect_perc_2008_crop)
east_aspect_perc_2010_vec<- as.vector(east_aspect_perc_2010_crop)
east_aspect_perc_2012_vec<- as.vector(east_aspect_perc_2012_crop)

south_aspect_perc_2008_vec<- as.vector(south_aspect_perc_2008_crop)
south_aspect_perc_2010_vec<- as.vector(south_aspect_perc_2010_crop)
south_aspect_perc_2012_vec<- as.vector(south_aspect_perc_2012_crop)

west_aspect_perc_2008_vec<- as.vector(west_aspect_perc_2008_crop)
west_aspect_perc_2010_vec<- as.vector(west_aspect_perc_2010_crop)
west_aspect_perc_2012_vec<- as.vector(west_aspect_perc_2012_crop)


## Calculate means and standard deviations

north_aspect_2008_mean<- mean(north_aspect_perc_2008_vec, na.rm = TRUE)
north_aspect_2008_sd<- sd(north_aspect_perc_2008_vec, na.rm = TRUE)
north_aspect_2010_mean<- mean(north_aspect_perc_2010_vec, na.rm = TRUE)
north_aspect_2010_sd<- sd(north_aspect_perc_2010_vec, na.rm = TRUE)
north_aspect_2012_mean<- mean(north_aspect_perc_2012_vec, na.rm = TRUE)
north_aspect_2012_sd<- sd(north_aspect_perc_2012_vec, na.rm = TRUE)

east_aspect_2008_mean<- mean(east_aspect_perc_2008_vec, na.rm = TRUE)
east_aspect_2008_sd<- sd(east_aspect_perc_2008_vec, na.rm = TRUE)
east_aspect_2010_mean<- mean(east_aspect_perc_2010_vec, na.rm = TRUE)
east_aspect_2010_sd<- sd(east_aspect_perc_2010_vec, na.rm = TRUE)
east_aspect_2012_mean<- mean(east_aspect_perc_2012_vec, na.rm = TRUE)
east_aspect_2012_sd<- sd(east_aspect_perc_2012_vec, na.rm = TRUE)

south_aspect_2008_mean<- mean(south_aspect_perc_2008_vec, na.rm = TRUE)
south_aspect_2008_sd<- sd(south_aspect_perc_2008_vec, na.rm = TRUE)
south_aspect_2010_mean<- mean(south_aspect_perc_2010_vec, na.rm = TRUE)
south_aspect_2010_sd<- sd(south_aspect_perc_2010_vec, na.rm = TRUE)
south_aspect_2012_mean<- mean(south_aspect_perc_2012_vec, na.rm = TRUE)
south_aspect_2012_sd<- sd(south_aspect_perc_2012_vec, na.rm = TRUE)

west_aspect_2008_mean<- mean(west_aspect_perc_2008_vec, na.rm = TRUE)
west_aspect_2008_sd<- sd(west_aspect_perc_2008_vec, na.rm = TRUE)
west_aspect_2010_mean<- mean(west_aspect_perc_2010_vec, na.rm = TRUE)
west_aspect_2010_sd<- sd(west_aspect_perc_2010_vec, na.rm = TRUE)
west_aspect_2012_mean<- mean(west_aspect_perc_2012_vec, na.rm = TRUE)
west_aspect_2012_sd<- sd(west_aspect_perc_2012_vec, na.rm = TRUE)



## Aspect bar plots with means and standard deviations

labels<- c("2008", "2010", "2012")
means<- c(north_aspect_2008_mean, north_aspect_2010_mean, north_aspect_2012_mean)
sds<- c(north_aspect_2008_sd, north_aspect_2010_sd, north_aspect_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - North Aspect", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 27, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(east_aspect_2008_mean, east_aspect_2010_mean, east_aspect_2012_mean)
sds<- c(east_aspect_2008_sd, east_aspect_2010_sd, east_aspect_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Percent Snow by East Aspect", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(south_aspect_2008_mean, south_aspect_2010_mean, south_aspect_2012_mean)
sds<- c(south_aspect_2008_sd, south_aspect_2010_sd, south_aspect_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - South Aspect", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(west_aspect_2008_mean, west_aspect_2010_mean, west_aspect_2012_mean)
sds<- c(west_aspect_2008_sd, west_aspect_2010_sd, west_aspect_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - West Aspect", col = c("#7570B3", "#1B9E77", "#D95F02"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


########### Analysis by vegetation class

east_veg_class<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/US_200EVT_ERarea/LF_fullarea/Lf_fullarea.tif')

class_dbf<- read.dbf('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/US_200EVT_ERarea/LF_fullarea/LF_fullarea.tif.vat.dbf')


## Call the raster attribute table ('RAT') for the raster data
east_veg_rat<- ratify(east_veg_class)

## Create a table that has 1 column and the same number of rows as levels of the data
rat<- levels(east_veg_rat)[[1]]
## Create a second column called 'landcover' with the names of the categories
rat$landcover<- c('Tree', 'Shrub', 'Herb', 'Water', 'Barren', 'Developed', 'Snow-Ice', 'Agriculture', 'Sparse')
## Assign the levels to the data table
levels(east_veg_rat)<- rat

## Change the projection of the shapefile so it can be displayed on top of the raster
veg_proj<- crs(east_veg_class)
ER_veg<- spTransform(ER, veg_proj)
proj4string(veg_proj)

## Plot the landcover classes (levels) of the data
levelplot(east_veg_rat, att = 'landcover', col.regions = c("#33A02C", "#B2DF8A", "#FDBF6F", "#1F78B4", "#B15928", "#6A3D9A","#A6CEE3", "#FFFF99", "#CAB2D6"), main = "Landcover classes") + layer(sp.polygons(ER_veg))

## Re-organize vegetation categories 
east_trees <- (east_veg_class == 1)
east_shrub <- (east_veg_class == 2)
east_grass <- (east_veg_class == 3 | east_veg_class == 8)  ## includes agriculture and herb areas
east_clear <- (east_veg_class == 5 | east_veg_class == 9) ## inlcudes sparse and barren areas
east_other<- (east_veg_class == 4 | east_veg_class == 6 | east_veg_class == 7) ## includes water, developed, and snow-ice areas


## Put the new categories into one raster
veg_class<- east_veg_class
veg_class[east_trees] <- 1
veg_class[east_shrub] <- 2
veg_class[east_grass] <- 3
veg_class[east_clear] <- 4
veg_class[east_other] <- 5


veg_class_fac<- as.factor(veg_class)
rat_veg<- levels(veg_class_fac)[[1]] 
rat_veg[["class"]] <- c("trees", "shrubs", "grass", "clear", "other")
levels(veg_class_fac)<- rat_veg

levelplot(veg_class_fac,  col.regions = c("#33A02C", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#6A3D9A"), main = "Landcover Classes") + layer(sp.polygons(ER_veg))


### Project the vegetation rasters to lat/long wgs84 (to match the dem/starfm data)
### Resample vegetation rasters to 30-m pixels (to match starfm data)
trees_proj<- projectRaster(east_trees, crs = data_proj, method = 'ngb')
trees_proj <- resample(trees_proj, snow_sum_2008, method = 'ngb')

shrubs_proj<- projectRaster(east_shrub, crs = data_proj, method = 'ngb')
shrubs_proj <- resample(shrubs_proj, snow_sum_2008, method = 'ngb')

grass_proj<- projectRaster(east_grass, crs = data_proj, method = 'ngb')
grass_proj <- resample(grass_proj, snow_sum_2008, method = 'ngb')

clear_proj<- projectRaster(east_clear, crs = data_proj, method = 'ngb')
clear_proj <- resample(clear_proj, snow_sum_2008, method = 'ngb')

other_proj<- projectRaster(east_other, crs = data_proj, method = 'ngb')
other_proj <- resample(other_proj, snow_sum_2008, method = 'ngb')

### Analyses by landcover type 

## Trees

trees_perc_2008 <- corrected_2008
masked <- Which(trees_proj == 0, cells = TRUE) ## identify all pixels that are not located in tree areas
trees_perc_2008[masked] <- NA ## mask out these pixels so only tree pixels have data

plot(trees_perc_2008, main = "2008 Annual Perc Snow for Tree Areas")
plot(ER_proj, border = 'black', add = TRUE)

trees_perc_2010 <- corrected_2010
masked <- Which(trees_proj == 0, cells = TRUE)
trees_perc_2010[masked] <- NA

plot(trees_perc_2010, main = "2010 Annual Perc Snow for Tree Areas")
plot(ER_proj, border = 'black', add = TRUE)

trees_perc_2012 <- corrected_2012
masked <- Which(trees_proj == 0, cells = TRUE)
trees_perc_2012[masked] <- NA

plot(trees_perc_2012, main = "2012 Annual Perc Snow for Tree Areas")
plot(ER_proj, border = 'black', add = TRUE)


## Shrubs

shrubs_perc_2008 <- corrected_2008
masked <- Which(shrubs_proj == 0, cells = TRUE) ## identify all pixels that are not located in shrub areas
shrubs_perc_2008[masked] <- NA ## mask out these pixels so only shrub pixels have data

plot(shrubs_perc_2008, main = "2008 Annual Perc Snow for Shrub Areas")
plot(ER_proj, border = 'black', add = TRUE)

shrubs_perc_2010 <- corrected_2010
masked <- Which(shrubs_proj == 0, cells = TRUE)
shrubs_perc_2010[masked] <- NA

plot(shrubs_perc_2010, main = "2010 Annual Perc Snow for Shrubs Areas")
plot(ER_proj, border = 'black', add = TRUE)

shrubs_perc_2012 <- corrected_2012
masked <- Which(shrubs_proj == 0, cells = TRUE)
shrubs_perc_2012[masked] <- NA

plot(shrubs_perc_2012, main = "2012 Annual Perc Snow for Shrubs Areas")
plot(ER_proj, border = 'black', add = TRUE)


## Grass

grass_perc_2008 <- corrected_2008
masked <- Which(grass_proj == 0, cells = TRUE) ## identify all pixels that are not located in grass areas
grass_perc_2008[masked] <- NA ## mask out these pixels so only grass pixels have data

plot(grass_perc_2008, main = "2008 Annual Perc Snow for Grass Areas")
plot(ER_proj, border = 'black', add = TRUE)

grass_perc_2010 <- corrected_2010
masked <- Which(grass_proj == 0, cells = TRUE)
grass_perc_2010[masked] <- NA

plot(grass_perc_2010, main = "2010 Annual Perc Snow for Grass Areas")
plot(ER_proj, border = 'black', add = TRUE)

grass_perc_2012 <- corrected_2012
masked <- Which(grass_proj == 0, cells = TRUE)
grass_perc_2012[masked] <- NA

plot(grass_perc_2012, main = "2012 Annual Perc Snow for Grass Areas")
plot(ER_proj, border = 'black', add = TRUE)


## Clear

clear_perc_2008 <- corrected_2008
masked <- Which(clear_proj == 0, cells = TRUE) ## identify all pixels that are not located in clear areas
clear_perc_2008[masked] <- NA ## mask out these pixels so only clear pixels have data
 
plot(clear_perc_2008, main = "2008 Annual Perc Snow for Clear Areas")
plot(ER_proj, border = 'black', add = TRUE)

clear_perc_2010 <- corrected_2010
masked <- Which(clear_proj == 0, cells = TRUE)
clear_perc_2010[masked] <- NA

plot(clear_perc_2010, main = "2010 Annual Perc Snow for Clear Areas")
plot(ER_proj, border = 'black', add = TRUE)

clear_perc_2012 <- corrected_2012
masked <- Which(clear_proj == 0, cells = TRUE)
clear_perc_2012[masked] <- NA

plot(clear_perc_2012, main = "2012 Annual Perc Snow for Clear Areas")
plot(ER_proj, border = 'black', add = TRUE)


## Other (did not use this category in analysis due to the small amount of data)

# other_perc_2008 <- corrected_2008
# masked <- Which(other_proj == 0, cells = TRUE)
# other_perc_2008[masked] <- NA
# 
# plot(other_perc_2008, main = "2008 Annual Perc Snow for Other Areas")
# plot(ER_proj, border = 'black', add = TRUE)
# 
# other_perc_2010 <- corrected_2010
# masked <- Which(other_proj == 0, cells = TRUE)
# other_perc_2010[masked] <- NA
# 
# plot(other_perc_2010, main = "2010 Annual Perc Snow for Other Areas")
# plot(ER_proj, border = 'black', add = TRUE)
# 
# other_perc_2012 <- corrected_2012
# masked <- Which(other_proj == 0, cells = TRUE)
# other_perc_2012[masked] <- NA
# 
# plot(other_perc_2012, main = "2012 Annual Perc Snow for Other Areas")
# plot(ER_proj, border = 'black', add = TRUE)

### Convert all rasters to vectors

trees_perc_2008_vec<- as.vector(trees_perc_2008)
trees_perc_2010_vec<- as.vector(trees_perc_2010)
trees_perc_2012_vec<- as.vector(trees_perc_2012)

shrubs_perc_2008_vec<- as.vector(shrubs_perc_2008)
shrubs_perc_2010_vec<- as.vector(shrubs_perc_2010)
shrubs_perc_2012_vec<- as.vector(shrubs_perc_2012)

grass_perc_2008_vec<- as.vector(grass_perc_2008)
grass_perc_2010_vec<- as.vector(grass_perc_2010)
grass_perc_2012_vec<- as.vector(grass_perc_2012)

clear_perc_2008_vec<- as.vector(clear_perc_2008)
clear_perc_2010_vec<- as.vector(clear_perc_2010)
clear_perc_2012_vec<- as.vector(clear_perc_2012)

# other_perc_2008_vec<- as.vector(other_perc_2008)
# other_perc_2010_vec<- as.vector(other_perc_2010)
# other_perc_2012_vec<- as.vector(other_perc_2012)


### Calclate means and standard deviations


trees_2008_mean<- mean(trees_perc_2008_vec, na.rm = TRUE)
trees_2008_sd<- sd(trees_perc_2008_vec, na.rm = TRUE)
trees_2010_mean<- mean(trees_perc_2010_vec, na.rm = TRUE)
trees_2010_sd<- sd(trees_perc_2010_vec, na.rm = TRUE)
trees_2012_mean<- mean(trees_perc_2012_vec, na.rm = TRUE)
trees_2012_sd<- sd(trees_perc_2012_vec, na.rm = TRUE)

shrubs_2008_mean<- mean(shrubs_perc_2008_vec, na.rm = TRUE)
shrubs_2008_sd<- sd(shrubs_perc_2008_vec, na.rm = TRUE)
shrubs_2010_mean<- mean(shrubs_perc_2010_vec, na.rm = TRUE)
shrubs_2010_sd<- sd(shrubs_perc_2010_vec, na.rm = TRUE)
shrubs_2012_mean<- mean(shrubs_perc_2012_vec, na.rm = TRUE)
shrubs_2012_sd<- sd(shrubs_perc_2012_vec, na.rm = TRUE)

grass_2008_mean<- mean(grass_perc_2008_vec, na.rm = TRUE)
grass_2008_sd<- sd(grass_perc_2008_vec, na.rm = TRUE)
grass_2010_mean<- mean(grass_perc_2010_vec, na.rm = TRUE)
grass_2010_sd<- sd(grass_perc_2010_vec, na.rm = TRUE)
grass_2012_mean<- mean(grass_perc_2012_vec, na.rm = TRUE)
grass_2012_sd<- sd(grass_perc_2012_vec, na.rm = TRUE)

clear_2008_mean<- mean(clear_perc_2008_vec, na.rm = TRUE)
clear_2008_sd<- sd(clear_perc_2008_vec, na.rm = TRUE)
clear_2010_mean<- mean(clear_perc_2010_vec, na.rm = TRUE)
clear_2010_sd<- sd(clear_perc_2010_vec, na.rm = TRUE)
clear_2012_mean<- mean(clear_perc_2012_vec, na.rm = TRUE)
clear_2012_sd<- sd(clear_perc_2012_vec, na.rm = TRUE)

# other_2008_mean<- mean(other_perc_2008_vec, na.rm = TRUE)
# other_2008_sd<- sd(other_perc_2008_vec, na.rm = TRUE)
# other_2010_mean<- mean(other_perc_2010_vec, na.rm = TRUE)
# other_2010_sd<- sd(other_perc_2010_vec, na.rm = TRUE)
# other_2012_mean<- mean(other_perc_2012_vec, na.rm = TRUE)
# other_2012_sd<- sd(other_perc_2012_vec, na.rm = TRUE)

## Landcover annual percent snow bar plot

labels<- c("2008", "2010", "2012")
means<- c(trees_2008_mean, trees_2010_mean, trees_2012_mean)
sds<- c(trees_2008_sd, trees_2010_sd, trees_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - Tree Areas", ylab="Mean Value", col = c("#7570B3", "#1B9E77", "#D95F02"), cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 24, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(shrubs_2008_mean, shrubs_2010_mean, shrubs_2012_mean)
sds<- c(shrubs_2008_sd, shrubs_2010_sd, shrubs_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - Shrub Areas", ylab="Mean Value",col = c("#7570B3", "#1B9E77", "#D95F02"), cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 22, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(grass_2008_mean, grass_2010_mean, grass_2012_mean)
sds<- c(grass_2008_sd, grass_2010_sd, grass_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - Grass Areas", ylab="Mean Value", col = c("#7570B3", "#1B9E77", "#D95F02"), cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 25, paste(round(means)), cex = 1.2)


labels<- c("2008", "2010", "2012")
means<- c(clear_2008_mean, clear_2010_mean, clear_2012_mean)
sds<- c(clear_2008_sd, clear_2010_sd, clear_2012_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Mean Annual Perc Snow - Clear Areas", ylab="Mean Value", col = c("#7570B3", "#1B9E77", "#D95F02"), cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3, las=1, ylim= range(pretty(c(0, 100))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 22, paste(round(means)), cex = 1.2)


# labels<- c("2008", "2010", "2012")
# means<- c(other_2008_mean, other_2010_mean, other_2012_mean)
# sds<- c(other_2008_sd, other_2010_sd, other_2012_sd)
# 
# par(mfrow=c(1,1))
# mids<- barplot(means, names.arg=labels, main="Mean Annual Percent Snow for Other Areas", ylab="Mean Value", col = brewer.pal(n = 3, name = "Dark2"), las=1, ylim= range(pretty(c(0, 100))))
# arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
# text(mids, means + 28, paste(round(means, digits = 2)), cex = 1)
