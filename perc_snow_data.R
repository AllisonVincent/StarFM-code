###### This script contains code for analyzing percent snow-covered area for a water year by landscape characteristics


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

snow_sum<- raster("./WY2012_snow_sum.tif") ## raster with the number of instances (days) each pixel was classified as "snow-covered" for the water year 2012

### Re-project the watershed shapefile to the projection of the data
data_proj<- crs(snow_sum)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)


## Plot the total number of snow covered days
dev.new(height=0.91*nrow(snow_sum)/50, width=1.09*ncol(snow_sum)/50)
par(mar = c(5,5,5,3.5))
plot(snow_sum, col = rev(topo.colors(140)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Total Snow Covered Days", cex = 2.0), line = 1.0)
plot(snow_sum, legend.only = TRUE, col = rev(topo.colors(140)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


######### Normalize the data by calculating the percent snow-covered days for the water year

total_model_sum<- raster('./WY2012_model_sum.tif') ## raster with the total number of instances (days) data were present for each pixel (regardless of snow status) for the water year 2012

corrected<- (snow_sum/total_model_sum) * 100 ## calculate the percent snow-covered days per-pixel for the water year


######### Compare number of snow-covered days to topographic features (Slope and aspect rasters created in dem_analysis.R script)

east_DEM<- raster('./East_DEM.tif')

east_slope<- raster('./east_slope.tif')

east_aspect<- raster('./east_aspect.tif')


######### Calculating elevation bands to use in analysis

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

low_elev_limit<- x[which(e_cdf >= 0.33)[1]] ## find the value of the 33rd percentile of data
mid_elev_limit<- x[which(e_cdf >= 0.66)[1]] ## find the value fo the 66th percentile of data


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

levelplot(elev_bands_fac, col.regions = terrain.colors(3), main = "Elevation Bands") + layer(sp.polygons(ER_proj))

### Save elevation band data to its own raster
#writeRaster(elev_bands, filename = "./elev_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)


############# Calculating slope bands to use in analysis

slope_stats<- stats(east_slope)
slope_max<- slope_stats[[8]]
slope_min<- slope_stats[[4]]


## Calculate the empirical cumulative distribution of the elevation data

slope_vec<- as.vector(east_slope)


x<- sort(slope_vec)

e_cdf<- 1:length(x)/ length(x)
#plot(x, e_cdf, type = 's')

low_slope_limit<- x[which(e_cdf >= 0.33)[1]] ## find the value of the 33rd percentile of data
mid_slope_limit<- x[which(e_cdf >= 0.66)[1]] ## find the value fo the 66th percentile of data


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

levelplot(slope_bands_fac, col.regions = rev(brewer.pal(3, name = "RdYlBu")), main = "Slope Bands") + layer(sp.polygons(ER_proj))

### Save slope band data to its own raster
#writeRaster(slope_bands, filename = "./slope_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)


########### Calculating aspect directions to use in analysis  


## Assign directions to variables

aspect_stats<- stats(east_aspect)

## Calculate aspect directions by degree value
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

levelplot(aspect_cat, att = 'direction', col.regions = c("#FF0000", "#FFA500", "#F0E68C", "#87CEEB"), main = 'Aspect') + layer(sp.polygons(ER_proj))


### Save aspect direction data to its own raster
#writeRaster(aspect_cat, filename = "./aspect_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)


########### Calculating land cover types to use in analysis 

# The raster of land cover type values
east_veg_class<- raster('./Lf_fullarea.tif')

# database file of above, needed for plotting the data categorically
class_dbf<- read.dbf('./LF_fullarea.tif.vat.dbf')


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


east_trees <- (east_veg_class == 1)
east_shrub <- (east_veg_class == 2)
east_grass <- (east_veg_class == 3 | east_veg_class == 8) ## includes agriculture and herb areas
east_clear <- (east_veg_class == 5 | east_veg_class == 9) ## includes sparse and barren areas
east_other<- (east_veg_class == 4 | east_veg_class == 6 | east_veg_class == 7) ## includes water, developed, and permanent snow-ice areas


## Same as above, just all in one raster
veg_class<- east_veg_class
veg_class[east_trees] <- 1
veg_class[east_shrub] <- 2
veg_class[east_grass] <- 3
veg_class[east_clear] <- 4
veg_class[east_other] <- 5


### Plot the new land cover types
veg_class_fac<- as.factor(veg_class)
rat_veg<- levels(veg_class_fac)[[1]] 
rat_veg[["class"]] <- c("trees", "shrubs", "grass", "clear", "other")
levels(veg_class_fac)<- rat_veg

levelplot(veg_class_fac,  col.regions = c("#33A02C", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#6A3D9A"), main = "Landcover Classes") + layer(sp.polygons(ER_veg))

### write new land cover type data to its own raster
#writeRaster(veg_class, filename = "./landcover_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)


### Project the landcover/vegetation rasters to lat/long wgs84 (to match starfm data)
### Resample landcover rasters to 30m pixels (to match starfm data)
trees_proj<- projectRaster(east_trees, crs = data_proj, method = 'ngb')
trees_proj <- resample(trees_proj, snow_sum, method = 'ngb')

shrubs_proj<- projectRaster(east_shrub, crs = data_proj, method = 'ngb')
shrubs_proj <- resample(shrubs_proj, snow_sum, method = 'ngb')

grass_proj<- projectRaster(east_grass, crs = data_proj, method = 'ngb')
grass_proj <- resample(grass_proj, snow_sum, method = 'ngb')

clear_proj<- projectRaster(east_clear, crs = data_proj, method = 'ngb')
clear_proj <- resample(clear_proj, snow_sum, method = 'ngb')

other_proj<- projectRaster(east_other, crs = data_proj, method = 'ngb')
other_proj <- resample(other_proj, snow_sum, method = 'ngb')



#################### Find the percentage of snow covered days by landscape attribute


##### By elevation

# For low elevations

perc_low_elev_snow <- corrected
masked <- Which(elev_low == 0, cells=TRUE) ## identify all pixels that are not located in the low elevation band
perc_low_elev_snow[masked] <- NA ## mask out those pixels so only low elevation areas have data


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_low_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Low Elev", cex = 1.5), line = 1.0)
plot(perc_low_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For mid elevations

perc_mid_elev_snow <- corrected
masked <- Which(elev_mid == 0, cells=TRUE) ## identify all pixels that are not located in the mid elevation band
perc_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elevation areas have data


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_mid_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Mid Elev", cex = 1.5), line = 1.0)
plot(perc_mid_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For high elevations

perc_high_elev_snow <- corrected
masked <- Which(elev_high == 0, cells=TRUE) ## identify all pixels that are not located in the high elevation band
perc_high_elev_snow[masked] <- NA ## mask out those pixels so only high elevation areas have data


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_high_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at High Elev", cex = 1.5), line = 1.0)
plot(perc_high_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


## box and whisker plot for percent snow covered days by elevation

low_elev_snow_vec<- as.vector(perc_low_elev_snow)

mid_elev_snow_vec<- as.vector(perc_mid_elev_snow)

high_elev_snow_vec<- as.vector(perc_high_elev_snow)


boxplot(low_elev_snow_vec, mid_elev_snow_vec, high_elev_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days at Elevation Band", cex.main = 1.2,
        at = c(1,2,3),
        names = c('Low', 'Mid', 'High'), cex.names = 1.2, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = "#D95F02",
        horizontal = TRUE)


######### By slope

# For low slopes

perc_low_slope_snow <- corrected
masked <- Which(slope_low == 0, cells=TRUE) ## identify all pixels that are not located in the low slope band
perc_low_slope_snow[masked] <- NA ## mask out those pixels so only low slope areas have data

## Need to crop raster due to edge effects from slope calculation
perc_low_slope_snow_crop<-crop(perc_low_slope_snow, extent(perc_low_slope_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_low_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Low Slope", cex = 1.5), line = 1.0)
plot(perc_low_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

# For mid slopes

perc_mid_slope_snow <- corrected
masked <- Which(slope_mid == 0, cells=TRUE) ## identify all pixels that are not located in the mid slope band
perc_mid_slope_snow[masked] <- NA ## mask out those pixels so only mid slope areas have data

## Need to crop raster due to edge effects from slope calcuation
perc_mid_slope_snow_crop<- crop(perc_mid_slope_snow, extent(perc_mid_slope_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_mid_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Mid Slope", cex = 1.5), line = 1.0)
plot(perc_mid_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

# For high slopes

perc_high_slope_snow <- corrected
masked <- Which(slope_high == 0, cells=TRUE) ## identify all pixels that are not located in the high slope band
perc_high_slope_snow[masked] <- NA ## mask out those pixels so only high slope areas have data

## Need to crop raster due to edge effects from slope calcuation
perc_high_slope_snow_crop<- crop(perc_high_slope_snow, extent(perc_high_slope_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_high_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at High Slope", cex = 1.5), line = 1.0)
plot(perc_high_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### Create box and whisker plot of percent snow cover by slope band

low_slope_snow_vec<- as.vector(perc_low_slope_snow_crop)

mid_slope_snow_vec<- as.vector(perc_mid_slope_snow_crop)

high_slope_snow_vec<- as.vector(perc_high_slope_snow_crop)


boxplot(low_slope_snow_vec, mid_slope_snow_vec, high_slope_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days at Slope Band", cex.main = 1.2,
        at = c(1,2,3),
        names = c('Low', 'Mid', 'High'), cex.names = 1.2, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = "#7570B3", 
        horizontal = TRUE)

############ By aspect


# For north facing slopes

perc_north_aspect_snow <- corrected
masked <- Which(north == 0, cells=TRUE) ## identify all pixels that are not north aspect
perc_north_aspect_snow[masked] <- NA ## mask out those pixels so only north aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_north_aspect_snow_crop<- crop(perc_north_aspect_snow, extent(perc_north_aspect_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_north_aspect_snow_crop, col = rev(topo.colors(150)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at North Aspect", cex = 1.5), line = 1.0)
plot(perc_north_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For east facing slopes

perc_east_aspect_snow <- corrected
masked <- Which(east == 0, cells=TRUE) ## identify all pixels that are not east aspect
perc_east_aspect_snow[masked] <- NA ## mask out those pixels so only east aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_east_aspect_snow_crop<- crop(perc_east_aspect_snow, extent(perc_east_aspect_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_east_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at East Aspect", cex = 1.5), line = 1.0)
plot(perc_east_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For south facing slopes

perc_south_aspect_snow <- corrected
masked <- Which(south == 0, cells=TRUE) ## identify all pixels that are not south aspect
perc_south_aspect_snow[masked] <- NA ## mask out those pixels so only south aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_south_aspect_snow_crop<- crop(perc_south_aspect_snow, extent(perc_south_aspect_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_south_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at South Aspect", cex = 1.5), line = 1.0)
plot(perc_south_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For west facing slopes

perc_west_aspect_snow <- corrected
masked <- Which(west == 0, cells=TRUE) ## identify all pixels that are not west aspect
perc_west_aspect_snow[masked] <- NA ## mask out those pixels so only west aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_west_aspect_snow_crop<- crop(perc_west_aspect_snow, extent(perc_west_aspect_snow, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_west_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at West Aspect", cex = 1.5), line = 1.0)
plot(perc_west_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### Create box and whisker plot of percent snow cover by aspect direction

north_aspect_snow_vec<- as.vector(perc_north_aspect_snow_crop)

east_aspect_snow_vec<- as.vector(perc_east_aspect_snow_crop)

south_aspect_snow_vec<- as.vector(perc_south_aspect_snow_crop)

west_aspect_snow_vec<- as.vector(perc_west_aspect_snow_crop)


boxplot(north_aspect_snow_vec, east_aspect_snow_vec, south_aspect_snow_vec, west_aspect_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days at Aspect", cex.main = 1.2,
        at = c(1,2,3,4),
        names = c('North', 'East', 'South', 'West'), cex.names = 1.1, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = ("#D95F02"), 
        horizontal = TRUE)


############## By landcover type


### For tree areas
perc_trees_snow <- corrected
masked <- Which(trees_proj == 0, cells = TRUE) ## identify all pixels that are not tree areas
perc_trees_snow[masked] <- NA ## mask out those pixels so only tree area pixels have data

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_trees_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Tree Areas", cex = 1.5), line = 1.0)
plot(perc_trees_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### For shrub areas
perc_shrubs_snow <- corrected
masked <- Which(shrubs_proj == 0, cells = TRUE) ## identify all pixels that are not shrub areas
perc_shrubs_snow[masked] <- NA ## mask out those pixels so only shrub area pixels have data

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_shrubs_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Shrub Areas", cex = 1.5), line = 1.0)
plot(perc_shrubs_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### For grass areas
perc_grass_snow <- corrected
masked <- Which(grass_proj == 0, cells = TRUE) ## identify all pixels that are not grass areas
perc_grass_snow[masked] <- NA ## mask out those pixels so only grass area pixels have data

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_grass_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Grass Areas", cex = 1.5), line = 1.0)
plot(perc_grass_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### For clear/unvegetated areas
perc_clear_snow <- corrected
masked <- Which(clear_proj == 0, cells=TRUE) ## identify all pixels that are not clear areas
perc_clear_snow[masked] <- NA ## mask out those pixels so only clear area pixels have data

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_clear_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Clear Areas", cex = 1.5), line = 1.0)
plot(perc_clear_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

### For all other areas
# perc_other_snow <- corrected
# masked <- Which(other_proj == 0, cells=TRUE)
# perc_other_snow[masked] <- NA
# 
# par(mar = c(2.5,2.5,2.5,2.5))
# plot(perc_other_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow for all other areas")
# plot(ER_proj, border = 'black', add = TRUE)


######### Create a box and whisker plot for percent snow by land cover type 

trees_snow_vec<- as.vector(perc_trees_snow)

shrubs_snow_vec<- as.vector(perc_shrubs_snow)

grass_snow_vec<- as.vector(perc_grass_snow)

clear_snow_vec<- as.vector(perc_clear_snow)

#other_snow_vec<- as.vector(perc_other_snow)


boxplot(trees_snow_vec, shrubs_snow_vec, grass_snow_vec, clear_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days by Landcover Class", cex.main = 1.2,
        at = c(1,2,3,4),
        names = c('Trees', 'Shrubs', 'Grass', 'Clear'), cex.names = 1.2, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = "#D95F02", 
        horizontal = TRUE)



############### Testing the hypothesis that elevation is the primary driver of snow cover patterns
##### Create box and whisker plots for slope by elevation band to isolate the elevation relationship for snow cover regardless of slope


## low slopes/low elevation
 
perc_low_slope_low_elev<- (elev_low == 1 & slope_low == 1) ## identify pixels that are located in both low elevation and low slope areas

perc_low_slope_low_elev_snow <- corrected
masked <- Which(perc_low_slope_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/low slope areas
perc_low_slope_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/low slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_low_slope_low_elev_snow_crop<- crop(perc_low_slope_low_elev_snow, extent(perc_low_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at low slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## low slopes/mid elevation

perc_low_slope_mid_elev<- (elev_mid == 1 & slope_low == 1) ## identify pixels that are located in both mid elevation and low slope areas

perc_low_slope_mid_elev_snow <- corrected
masked <- Which(perc_low_slope_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/low slope areas
perc_low_slope_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/low slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_low_slope_mid_elev_snow_crop<- crop(perc_low_slope_mid_elev_snow, extent(perc_low_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at low slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## low slopes/high elevation

perc_low_slope_high_elev<- (elev_high == 1 & slope_low == 1) ## identify pixels that are located in both high elevation and low slope areas

perc_low_slope_high_elev_snow <- corrected
masked <- Which(perc_low_slope_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/low slope areas
perc_low_slope_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/low slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_low_slope_high_elev_snow_crop<- crop(perc_low_slope_high_elev_snow, extent(perc_low_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at low slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    


### create a box and whisker plot for low slope and all elevation relationships

perc_low_slope_low_elev_snow_crop_vec<- as.vector(perc_low_slope_low_elev_snow_crop)

perc_low_slope_mid_elev_snow_crop_vec<- as.vector(perc_low_slope_mid_elev_snow_crop)

perc_low_slope_high_elev_snow_crop_vec<- as.vector(perc_low_slope_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_low_slope_low_elev_snow_crop_vec, perc_low_slope_mid_elev_snow_crop_vec, perc_low_slope_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered Days at Low Slopes/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#91BFDB", 
        horizontal = TRUE)



## mid slopes/low elevation

perc_mid_slope_low_elev<- (elev_low == 1 & slope_mid == 1) ## identify pixels that are located in both low elevation and mid slope areas

perc_mid_slope_low_elev_snow <- corrected
masked <- Which(perc_mid_slope_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/mid slope areas
perc_mid_slope_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/mid slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_mid_slope_low_elev_snow_crop<- crop(perc_mid_slope_low_elev_snow, extent(perc_mid_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## mid slopes/mid elevation

perc_mid_slope_mid_elev<- (elev_mid == 1 & slope_mid == 1) ## identify pixels that are located in both mid elevation and mid slope areas

perc_mid_slope_mid_elev_snow <- corrected
masked <- Which(perc_mid_slope_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/mid slope areas
perc_mid_slope_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/mid slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_mid_slope_mid_elev_snow_crop<- crop(perc_mid_slope_mid_elev_snow, extent(perc_mid_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## mid slopes/high elevation

perc_mid_slope_high_elev<- (elev_high == 1 & slope_mid == 1) ## identify pixels that are located in both high elevation and mid slope areas

perc_mid_slope_high_elev_snow <- corrected
masked <- Which(perc_mid_slope_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/mid slope areas
perc_mid_slope_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/mid slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_mid_slope_high_elev_snow_crop<- crop(perc_mid_slope_high_elev_snow, extent(perc_mid_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    



### create a box and whisker plot for mid slope and all elevation relationships

perc_mid_slope_low_elev_snow_crop_vec<- as.vector(perc_mid_slope_low_elev_snow_crop)

perc_mid_slope_mid_elev_snow_crop_vec<- as.vector(perc_mid_slope_mid_elev_snow_crop)

perc_mid_slope_high_elev_snow_crop_vec<- as.vector(perc_mid_slope_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_mid_slope_low_elev_snow_crop_vec, perc_mid_slope_mid_elev_snow_crop_vec, perc_mid_slope_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered Days at Mid Slopes/All Elevations", cex.main  = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), ces.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#FFFFBF", 
        horizontal = TRUE)



## high slopes/low elevation

perc_high_slope_low_elev<- (elev_low == 1 & slope_high == 1) ## identify pixels that are located in both low elevation and high slope areas

perc_high_slope_low_elev_snow <- corrected
masked <- Which(perc_high_slope_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/high slope areas
perc_high_slope_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/high slope pixels have data

# Need to crop raster due to edge effects from slope calcuation
perc_high_slope_low_elev_snow_crop<- crop(perc_high_slope_low_elev_snow, extent(perc_high_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## high slopes/mid elevation

perc_high_slope_mid_elev<- (elev_mid == 1 & slope_high == 1) ## identify pixels that are located in both mid elevation and high slope areas

perc_high_slope_mid_elev_snow <- corrected
masked <- Which(perc_high_slope_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/high slope areas
perc_high_slope_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/high slope pixels have data

## Need to crop raster due to edge effects from slope calcuation
perc_high_slope_mid_elev_snow_crop<- crop(perc_high_slope_mid_elev_snow, extent(perc_high_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## high slopes/high elevation

perc_high_slope_high_elev<- (elev_high == 1 & slope_high == 1) ## identify pixels that are located in both high elevation and high slope areas

perc_high_slope_high_elev_snow <- corrected
masked <- Which(perc_high_slope_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/high slope areas
perc_high_slope_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/high slope pixels have data

## Need to crop raster due to edge effects from slope calcuation
perc_high_slope_high_elev_snow_crop<- crop(perc_high_slope_high_elev_snow, extent(perc_high_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    


### create a box and whisker plot for high slope and all elevation relationships

perc_high_slope_low_elev_snow_crop_vec<- as.vector(perc_high_slope_low_elev_snow_crop)

perc_high_slope_mid_elev_snow_crop_vec<- as.vector(perc_high_slope_mid_elev_snow_crop)

perc_high_slope_high_elev_snow_crop_vec<- as.vector(perc_high_slope_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_high_slope_low_elev_snow_crop_vec, perc_high_slope_mid_elev_snow_crop_vec, perc_high_slope_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered days at High Slopes/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#FC8D59", 
        horizontal = TRUE)


###### Create boxplots for aspect by elevation band to isolate the elevation relationship between snow cover and landscape characteristics

## north aspect/low elevation

perc_north_aspect_low_elev<- (elev_low == 1 & north == 1) ## identify pixels that are located in both low elevation and north aspect areas

perc_north_aspect_low_elev_snow <- corrected
masked <- Which(perc_north_aspect_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/north aspect areas
perc_north_aspect_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/north aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_north_aspect_low_elev_snow_crop<- crop(perc_north_aspect_low_elev_snow, extent(perc_north_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## north aspect/mid elevation

perc_north_aspect_mid_elev<- (elev_mid == 1 & north == 1) ## identify pixels that are located in both mid elevation and north aspect areas

perc_north_aspect_mid_elev_snow <- corrected
masked <- Which(perc_north_aspect_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/north aspect areas
perc_north_aspect_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/north aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_north_aspect_mid_elev_snow_crop<- crop(perc_north_aspect_mid_elev_snow, extent(perc_north_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## north aspect/high elevation

perc_north_aspect_high_elev<- (elev_high == 1 & north == 1) ## identify pixels that are located in both high elevation and north aspect areas

perc_north_aspect_high_elev_snow <- corrected
masked <- Which(perc_north_aspect_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/north aspect areas
perc_north_aspect_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/north aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_north_aspect_high_elev_snow_crop<- crop(perc_north_aspect_high_elev_snow, extent(perc_north_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for north aspect and all elevation relationships

perc_north_aspect_low_elev_snow_crop_vec<- as.vector(perc_north_aspect_low_elev_snow_crop)

perc_north_aspect_mid_elev_snow_crop_vec<- as.vector(perc_north_aspect_mid_elev_snow_crop)

perc_north_aspect_high_elev_snow_crop_vec<- as.vector(perc_north_aspect_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_north_aspect_low_elev_snow_crop_vec, perc_north_aspect_mid_elev_snow_crop_vec, perc_north_aspect_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered Days at North Aspect/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#FF0000", 
        horizontal = TRUE)



## east aspect/low elevation

perc_east_aspect_low_elev<- (elev_low == 1 & east == 1) ## identify pixels that are located in both low elevation and east aspect areas

perc_east_aspect_low_elev_snow <- corrected
masked <- Which(perc_east_aspect_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/east aspect areas
perc_east_aspect_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/east aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_east_aspect_low_elev_snow_crop<- crop(perc_east_aspect_low_elev_snow, extent(perc_east_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## east aspect/mid elevation

perc_east_aspect_mid_elev<- (elev_mid == 1 & east == 1) ## identify pixels that are located in both mid elevation and east aspect areas

perc_east_aspect_mid_elev_snow <- corrected
masked <- Which(perc_east_aspect_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/east aspect areas
perc_east_aspect_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/east aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_east_aspect_mid_elev_snow_crop<- crop(perc_east_aspect_mid_elev_snow, extent(perc_east_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## east aspect/high elevation

perc_east_aspect_high_elev<- (elev_high == 1 & east == 1) ## identify pixels that are located in both high elevation and east aspect areas

perc_east_aspect_high_elev_snow <- corrected
masked <- Which(perc_east_aspect_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/east aspect areas
perc_east_aspect_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/east aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_east_aspect_high_elev_snow_crop<- crop(perc_east_aspect_high_elev_snow, extent(perc_east_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for east aspect and all elevation relationships

perc_east_aspect_low_elev_snow_crop_vec<- as.vector(perc_east_aspect_low_elev_snow_crop)

perc_east_aspect_mid_elev_snow_crop_vec<- as.vector(perc_east_aspect_mid_elev_snow_crop)

perc_east_aspect_high_elev_snow_crop_vec<- as.vector(perc_east_aspect_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_east_aspect_low_elev_snow_crop_vec, perc_east_aspect_mid_elev_snow_crop_vec, perc_east_aspect_high_elev_snow_crop_vec, 
        main = "WY 2008 Percent Snow-Covered Days at East Aspect/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,  
        col = "#FFA500", 
        horizontal = TRUE)


## south aspect/low elevation

perc_south_aspect_low_elev<- (elev_low == 1 & south == 1) ## identify pixels that are located in both low elevation and south aspect areas

perc_south_aspect_low_elev_snow <- corrected
masked <- Which(perc_south_aspect_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/south aspect areas
perc_south_aspect_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/south aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_south_aspect_low_elev_snow_crop<- crop(perc_south_aspect_low_elev_snow, extent(perc_south_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## south aspect/mid elevation

perc_south_aspect_mid_elev<- (elev_mid == 1 & south == 1) ## identify pixels that are located in both mid elevation and south aspect areas

perc_south_aspect_mid_elev_snow <- corrected
masked <- Which(perc_south_aspect_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/south aspect areas
perc_south_aspect_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/south aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_south_aspect_mid_elev_snow_crop<- crop(perc_south_aspect_mid_elev_snow, extent(perc_south_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## south aspect/high elevation

perc_south_aspect_high_elev<- (elev_high == 1 & south == 1) ## identify pixels that are located in both high elevation and south aspect areas

perc_south_aspect_high_elev_snow <- corrected
masked <- Which(perc_south_aspect_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/south aspect areas
perc_south_aspect_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/south aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_south_aspect_high_elev_snow_crop<- crop(perc_south_aspect_high_elev_snow, extent(perc_south_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for south aspect and all elevation relationships

perc_south_aspect_low_elev_snow_crop_vec<- as.vector(perc_south_aspect_low_elev_snow_crop)

perc_south_aspect_mid_elev_snow_crop_vec<- as.vector(perc_south_aspect_mid_elev_snow_crop)

perc_south_aspect_high_elev_snow_crop_vec<- as.vector(perc_south_aspect_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_south_aspect_low_elev_snow_crop_vec, perc_south_aspect_mid_elev_snow_crop_vec, perc_south_aspect_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered Days at South Aspect/All Elevations", cex.main = 1.0, 
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2, 
        col = "#F0E68C", 
        horizontal = TRUE)


## west aspect/low elevation

perc_west_aspect_low_elev<- (elev_low == 1 & west == 1) ## identify pixels that are located in both low elevation and west aspect areas

perc_west_aspect_low_elev_snow <- corrected
masked <- Which(perc_west_aspect_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/west aspect areas
perc_west_aspect_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/west aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_west_aspect_low_elev_snow_crop<- crop(perc_west_aspect_low_elev_snow, extent(perc_west_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## west aspect/mid elevation

perc_west_aspect_mid_elev<- (elev_mid == 1 & west == 1) ## identify pixels that are located in both mid elevation and west aspect areas

perc_west_aspect_mid_elev_snow <- corrected
masked <- Which(perc_west_aspect_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/west aspect areas
perc_west_aspect_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/west aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_west_aspect_mid_elev_snow_crop<- crop(perc_west_aspect_mid_elev_snow, extent(perc_west_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## west aspect/high elevation

perc_west_aspect_high_elev<- (elev_high == 1 & west == 1) ## identify pixels that are located in both high elevation and west aspect areas

perc_west_aspect_high_elev_snow <- corrected
masked <- Which(perc_west_aspect_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/west aspect areas
perc_west_aspect_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/west aspect pixels have data

## Need to crop raster due to edge effects from aspect calcuation
perc_west_aspect_high_elev_snow_crop<- crop(perc_west_aspect_high_elev_snow, extent(perc_west_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for west aspect and all elevation relationships

perc_west_aspect_low_elev_snow_crop_vec<- as.vector(perc_west_aspect_low_elev_snow_crop)

perc_west_aspect_mid_elev_snow_crop_vec<- as.vector(perc_west_aspect_mid_elev_snow_crop)

perc_west_aspect_high_elev_snow_crop_vec<- as.vector(perc_west_aspect_high_elev_snow_crop)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_west_aspect_low_elev_snow_crop_vec, perc_west_aspect_mid_elev_snow_crop_vec, perc_west_aspect_high_elev_snow_crop_vec, 
        main = "WY 2012 Percent Snow-Covered Days at West Aspect/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0, 
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#87CEEB", 
        horizontal = TRUE)


###### Create boxplots for land cover type by elevation band to isolate the elevation relationship between snow cover and landscape characteristics


## tree areas/low elevation

perc_trees_low_elev<- (elev_low == 1 & trees_proj == 1) ## identify pixels that are located in both low elevation and tree areas

perc_trees_low_elev_snow <- corrected
masked <- Which(perc_trees_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/tree areas
perc_trees_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/tree pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## tree areas/mid elevation

perc_trees_mid_elev<- (elev_mid == 1 & trees_proj == 1) ## identify pixels that are located in both mid elevation and tree areas

perc_trees_mid_elev_snow <- corrected
masked <- Which(perc_trees_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/tree areas
perc_trees_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/tree pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## tree areas/high elevation

perc_trees_high_elev<- (elev_high == 1 & trees_proj == 1) ## identify pixels that are located in both high elevation and tree areas

perc_trees_high_elev_snow <- corrected
masked <- Which(perc_trees_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/tree areas
perc_trees_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/tree pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for tree areas and all elevation relationships

perc_trees_low_elev_snow_vec<- as.vector(perc_trees_low_elev_snow)

perc_trees_mid_elev_snow_vec<- as.vector(perc_trees_mid_elev_snow)

perc_trees_high_elev_snow_vec<- as.vector(perc_trees_high_elev_snow)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_trees_low_elev_snow_vec, perc_trees_mid_elev_snow_vec, perc_trees_high_elev_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days for Tree Areas/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#33A02C", 
        horizontal = TRUE)


## shrub areas/low elevation

perc_shrubs_low_elev<- (elev_low == 1 & shrubs_proj == 1) ## identify pixels that are located in both low elevation and shrub areas

perc_shrubs_low_elev_snow <- corrected
masked <- Which(perc_shrubs_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/shrub areas
perc_shrubs_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/shrub pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## shrub areas/mid elevation

perc_shrubs_mid_elev<- (elev_mid == 1 & shrubs_proj == 1) ## identify pixels that are located in both mid elevation and shrub areas

perc_shrubs_mid_elev_snow <- corrected
masked <- Which(perc_shrubs_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/shrub areas
perc_shrubs_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/shrub pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## shrub areas/high elevation

perc_shrubs_high_elev<- (elev_high == 1 & shrubs_proj == 1) ## identify pixels that are located in both high elevation and shrub areas

perc_shrubs_high_elev_snow <- corrected
masked <- Which(perc_shrubs_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/shrub areas
perc_shrubs_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/shrub pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for shrub areas and all elevation relationships

perc_shrubs_low_elev_snow_vec<- as.vector(perc_shrubs_low_elev_snow)

perc_shrubs_mid_elev_snow_vec<- as.vector(perc_shrubs_mid_elev_snow)

perc_shrubs_high_elev_snow_vec<- as.vector(perc_shrubs_high_elev_snow)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_shrubs_low_elev_snow_vec, perc_shrubs_mid_elev_snow_vec, perc_shrubs_high_elev_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days for Shrub Areas/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#B2DF8A", 
        horizontal = TRUE)


## grass areas/low elevation

perc_grass_low_elev<- (elev_low == 1 & grass_proj == 1) ## identify pixels that are located in both low elevation and grass areas

perc_grass_low_elev_snow <- corrected
masked <- Which(perc_grass_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/grass areas
perc_grass_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/grass pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## grass areas/mid elevation

perc_grass_mid_elev<- (elev_mid == 1 & grass_proj == 1) ## identify pixels that are located in both mid elevation and grass areas

perc_grass_mid_elev_snow <- corrected
masked <- Which(perc_grass_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/grass areas
perc_grass_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/grass pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## grass areas/high elevation

perc_grass_high_elev<- (elev_high == 1 & grass_proj == 1) ## identify pixels that are located in both high elevation and grass areas

perc_grass_high_elev_snow <- corrected
masked <- Which(perc_grass_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/grass areas
perc_grass_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/grass pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for grass areas and all elevation relationships

perc_grass_low_elev_snow_vec<- as.vector(perc_grass_low_elev_snow)

perc_grass_mid_elev_snow_vec<- as.vector(perc_grass_mid_elev_snow)

perc_grass_high_elev_snow_vec<- as.vector(perc_grass_high_elev_snow)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_grass_low_elev_snow_vec, perc_grass_mid_elev_snow_vec, perc_grass_high_elev_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days for Grass Areas/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#CAB2D6", 
        horizontal = TRUE)


## clear areas/low elevation

perc_clear_low_elev<- (elev_low == 1 & clear_proj == 1) ## identify pixels that are located in both low elevation and clear areas

perc_clear_low_elev_snow <- corrected
masked <- Which(perc_clear_low_elev == 0, cells=TRUE) ## identify all pixels that are not low elev/clear areas
perc_clear_low_elev_snow[masked] <- NA ## mask out those pixels so only low elev/clear pixels have data

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## clear areas/mid elevation

perc_clear_mid_elev<- (elev_mid == 1 & clear_proj == 1) ## identify pixels that are located in both mid elevation and clear areas

perc_clear_mid_elev_snow <- corrected
masked <- Which(perc_clear_mid_elev == 0, cells=TRUE) ## identify all pixels that are not mid elev/clear areas
perc_clear_mid_elev_snow[masked] <- NA ## mask out those pixels so only mid elev/clear pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## clear areas/high elevation

perc_clear_high_elev<- (elev_high == 1 & clear_proj == 1) ## identify pixels that are located in both high elevation and clear areas

perc_clear_high_elev_snow <- corrected
masked <- Which(perc_clear_high_elev == 0, cells=TRUE) ## identify all pixels that are not high elev/clear areas
perc_clear_high_elev_snow[masked] <- NA ## mask out those pixels so only high elev/clear pixels have data


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create a box and whisker plot for clear areas and all elevation relationships

perc_clear_low_elev_snow_vec<- as.vector(perc_clear_low_elev_snow)

perc_clear_mid_elev_snow_vec<- as.vector(perc_clear_mid_elev_snow)

perc_clear_high_elev_snow_vec<- as.vector(perc_clear_high_elev_snow)


par(mar = c(4.5,4.5,4.5,3.5))
boxplot(perc_clear_low_elev_snow_vec, perc_clear_mid_elev_snow_vec, perc_clear_high_elev_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days for Clear Areas/All Elevations", cex.main = 1.0,
        at = c(1,2,3),
        names = c('Low elev', 'Mid elev', 'High elev'), cex.names = 1.2, cex.axis = 1.0,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.2,
        col = "#6A3D9A", 
        horizontal = TRUE)

