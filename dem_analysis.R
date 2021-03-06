## This script is for calculating landscape characteristics of slope and aspect from DEM data

## DEM data for this project is STRM 30m data downloaded from Google Earth Engine

library(sp)
library(sf)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(caret)
library(data.table)
library(wesanderson)

### Load the DEM and watershed shapefile, set to common projection

east_dem <- raster('./East_DEM.tif')
ER<- readOGR('./EastRiver_Project.shp') #shapefile of the study watershed for reference

## Find the projection of the DEM and project the watershed shapefile to match
proj4string(east_dem)
data_proj<- crs(east_dem)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

### Plot the data 

plot(east_dem, main = "Elevation (m)")
plot(ER_proj, border = "black", add = TRUE)


### Find slope direction and aspect
## Slope computed using the Horn (1981) method 
east_slope<- terrain(east_dem, opt = "slope", unit = "degrees", neighbors = 8)

plot(east_slope, main = "Slope (deg)")
plot(ER_proj, border = "black", add = TRUE)

## Option to save slope results as its own raster file
#writeRaster(east_slope, filename="./east_slope.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)

east_aspect<- terrain(east_dem, opt = "aspect", unit = "degrees", neighbors = 8)

plot(east_aspect, main = "Aspect (deg)")
plot(ER_proj, border = "black", add = TRUE)

## Option to save aspect results as its own raster file
#writeRaster(east_aspect, filename="./east_aspect.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)


###################################################################################

### This section of code is for calculating elevation statistics for only the area within the watershed boundary

## crop raster using watershed vector

## need to use the mask function, otherwise crop will keep the area a square of the largest extents of the watershed boundary
east_dem_mask<- mask(east_dem, ER_proj, inverse = FALSE)
plot(east_dem_mask)

mask_stats<- stats(east_dem_mask)
mask_min<- mask_stats[[4]]
mask_max<- mask_stats[[8]]

mask_range<- mask_max-mask_min

## Find the cumulative area in km^2 for just the watershed
area<- cellStats(area(east_dem_mask, na.rm = TRUE), stat = sum)
