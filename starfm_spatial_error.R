###### This script can be used to identify STARFM performance spatially. Here confusion matrix metrics (true positives, true negatives, false positives, false negatives) are computed for Landsat and STARFM data manually to preserve the spatial information of where these metrics occur. 

#### Inputs for this script are generated from the 


setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season')


library(sp)
library(sf)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(caret)
library(data.table)
library(fields)


landsat<- raster("./10.16.2015/landsat_testday_101615.tif")
starfm<- raster("./10.16.2015/starfm_testday_101615.tif")
ER<- readOGR("./ER_watershed_shp/EastRiver_Project.shp") #shapefile of the study watershed for reference

data_proj<- crs(landsat)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

#landsat<- brick(landsat_all[[10]])
#starfm<- brick(starfm_all[[10]])


## Isolating and plotting the NDSI for the test dates

#landsat_test<- brick(landsat[[8]])
landsat_masked<- Which(landsat == 0, cells = TRUE)
landsat[landsat_masked] <- -11111
landsat_plot<- reclassify(landsat,cbind(-Inf, -11111, NA))

#starfm_test<- brick(starfm[[17]])
starfm_plot <- reclassify(starfm,cbind(-Inf, -11111, NA))

par(mfrow=c(1,2))
plot(landsat_plot, col = rev(cm.colors(40)), main = "NDSI from Landsat")
plot(ER_proj, border = "black", add = TRUE)
plot(starfm_plot, col = rev(cm.colors(40)), main = "NDSI from Model")
plot(ER_proj, border = "black", add = TRUE)

#writeRaster(landsat, filename="./12.28.2015/landsat_testday_122815.tif", bandorder='BSQ', 
#          datatype='INT2S',format='GTiff', overwrite=TRUE)


#writeRaster(starfm, filename="./12.28.2015/starfm_testday_122815.tif", bandorder='BSQ', 
#            datatype='INT2S',format='GTiff', overwrite=TRUE)




## Creating binary snow rasters

# Any values in our rasters of below 4000 are given a new value of zero
landsat_binary<- reclassify(landsat_plot,cbind(-Inf, 4000, 0))
starfm_binary<- reclassify(starfm_plot,cbind(-Inf, 4000, 0))

# Any values in our rasters of 4000 or above are given a new value of 1
landsat_binary<- reclassify(landsat_binary,cbind(4000, Inf, 1))
starfm_binary<- reclassify(starfm_binary,cbind(4000, Inf, 1))

par(mfrow=c(1,2)) 
plot(landsat_binary,  col = rev(cm.colors(40)), main = "Snow/No-Snow from Landsat")
plot(ER_proj, border = "black", add = TRUE)
plot(starfm_binary, col = rev(cm.colors(40)), main = "Snow/No-Snow from Model")
plot(ER_proj, border = "black", add = TRUE)

## Manually compare snow/no-snow between the two rasters to preserve the spatial location of each pixel

tp<- (landsat_binary == 1 & starfm_binary == 1) 
tn<- (landsat_binary == 0 & starfm_binary == 0)
fp<- (landsat_binary == 0 & starfm_binary == 1)
fn<- (landsat_binary == 1 & starfm_binary == 0)

## Plot to show spatially where the correct predictions and errors are

par(mfrow=c(2,2))
plot(tp, main = "True Positives 10.16.2015")
plot(ER_proj, border = "black", add = TRUE)
plot(tn, main = "True Negatives")
plot(ER_proj, border = "black", add = TRUE)
plot(fp, main = "False Positives")
plot(ER_proj, border = "black", add = TRUE)
plot(fn, main = "False Negatives")
plot(ER_proj, border = "black", add = TRUE)

## Make plots showing all areas that were predicted correctly and those that were not

good<- (tp == 1 | tn == 1)
bad<- (fp == 1 | fn == 1)

par(mfrow=c(1,2))
plot(good, col=topo.colors(3), main = "Correctly Identified")
plot(ER_proj, border = "black", add = TRUE)
plot(bad, col=rev(topo.colors(3)), main = "Misidentified")
plot(ER_proj, border = "black", add = TRUE)

par(mfrow=c(1,1))
plot(good, col=topo.colors(3), main = "Correctly Identified 10.16.2015")
plot(ER_proj, border = "black", add = TRUE)


#writeRaster(good, filename="./12.28.2015/id_122815.tif", bandorder='BSQ', 
#            datatype='INT2S',format='GTiff', overwrite=TRUE)


#writeRaster(good, filename="./model_id/id_122815.tif", bandorder='BSQ', 
#            datatype='INT2S',format='GTiff', overwrite=TRUE)


## Create a raster brick that has the number of instances there was a model prediction at each pixel 

## Put all the rasters of the id'd dates (rasters of 1 for correct prediction, 0 for incorrect prediction), into one brick raster

setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/model_id')

brick_id<- do.call(brick, lapply(list.files(path = getwd(), pattern = "id_*"), raster))

## Convert all pixels that have a prediction to 1 so they can be summed, keep all pixels that didn't have a prediction as NA

brick_ones <- brick_id

brick_ones[brick_ones == 0] <- 1

## Sum the total number of model predictions per pixel (regardless of correctness)

model_pred_sum<- calc(brick_ones, fun = sum, na.rm = TRUE)


plot(model_pred_sum, main = "Sums of all model predictions per pixels")
plot(ER_proj, border = "black", add = TRUE)


#writeRaster(model_pred_sum, filename="./model_id/model_pred_sum.tif", bandorder='BSQ', 
#            datatype='INT2S',format='GTiff', overwrite=TRUE)
