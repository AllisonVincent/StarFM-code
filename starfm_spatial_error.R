###### This script can be used to identify STARFM performance spatially. Here confusion matrix metrics (true positives, true negatives, false positives, false negatives) are computed for Landsat and STARFM data manually to preserve the spatial information of where these metrics occur. 

#### Inputs for this script (binary snow status rasters) are generated from the conf_matrix.R script


library(sp)
library(sf)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(caret)
library(data.table)
library(fields)


landsat_binary<- raster("./landsat_binary.tif")
starfm_binary<- raster("./starfm_binary.tif")
ER<- readOGR("./ER_watershed_shp/EastRiver_Project.shp") #shapefile of the study watershed for reference

## Project the watershed shapefile to match the projection of the data
data_proj<- crs(landsat_binary)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)


##################### Manually compare snow/no-snow between the two rasters to preserve the spatial location of each pixel

tp<- (landsat_binary == 1 & starfm_binary == 1) 
tn<- (landsat_binary == 0 & starfm_binary == 0)
fp<- (landsat_binary == 0 & starfm_binary == 1)
fn<- (landsat_binary == 1 & starfm_binary == 0)

##### Plot to show spatially where the correct predictions and errors are
## A value of 1 in the plots below means the metric of interest is present, a value of 0 means that it is not

par(mfrow=c(2,2))
plot(tp, main = "True Positives")
plot(ER_proj, border = "black", add = TRUE) ## the locations of all the TPs (values of 1)
plot(tn, main = "True Negatives")
plot(ER_proj, border = "black", add = TRUE) ## the locations of all the TNs (values of 1)
plot(fp, main = "False Positives")
plot(ER_proj, border = "black", add = TRUE) ## the locations of all the FPs (values of 1)
plot(fn, main = "False Negatives")
plot(ER_proj, border = "black", add = TRUE) ## the locations of all the FNs (values of 1)

## Make plots showing all areas that were predicted correctly and those that were not

good<- (tp == 1 | tn == 1) ## all pixels correctly predicted
bad<- (fp == 1 | fn == 1) ## all pixels incorrectly predicted (can be plotted for reference)

par(mfrow=c(1,2))
plot(good, col=rev(topo.colors(3)), main = "Correctly Identified") ## value of 1 = correctly identified
plot(ER_proj, border = "black", add = TRUE)


### Save the "good" raster file with the results of the date being analyzed. Values of 1 in this raster equal a correct prediction, values of 0 equal an incorrect prediction
writeRaster(good, filename="./id_date_of_observation.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)


############## Create a raster brick that has the number of correct predictions per pixel, with one date of analysis per layer

## Save all rasters generated above that have 'id' in the front of their name in the same folder, then can create a single raster brick out of them all into one brick raster

## setwd()  ## call the folder with the rasters here

## combine all rasters with "id" prefix
brick_id<- do.call(brick, lapply(list.files(path = getwd(), pattern = "id_*"), raster))

## create a raster brick from the combined raster files
all_sum<- calc(brick, fun = sum, na.rm = TRUE)

writeRaster(all_sum, filename="./east_sum.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)

## Create another raster brick that has the number of instances there was a model prediction at each pixel (regardless of correctness), with one date of analysis per layer

## Convert all pixels that have a prediction to 1 so they can be summed, keep all pixels that didn't have a prediction as NA

brick_ones <- brick_id

brick_ones[brick_ones == 0] <- 1

## Sum the total number of model predictions per pixel (regardless of correctness)

model_pred_sum<- calc(brick_ones, fun = sum, na.rm = TRUE)


plot(model_pred_sum, main = "Sums of all model predictions per pixels")
plot(ER_proj, border = "black", add = TRUE)


writeRaster(model_pred_sum, filename="./model_pred_sum.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)
