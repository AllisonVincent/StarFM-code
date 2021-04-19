#### This script is for the purpose of viewing individual rasters of STARM and Landsat NDSI data, as well as creating a confusion matrix to measure the performance of STARFM against the Landsat for model validation.


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
library(SpatialTools)
library(foreign)
library(oceanmap)
library(tidyr)
library(plotly)
library(rasterVis)


########## Isolating and plotting the NDSI for the test dates

## Load the data files needed
landsat<- brick("./landsat.tif")
modis<- brick("./mod.tif")
starfm<- brick("./starfm.tif")
ER<- readOGR("./EastRiver_Project.shp") #shapefile of the study watershed for reference


## Project the shapefile to match the raster files 
data_proj<- crs(starfm)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

###### This section contains code to isolate the raster layers for a specific date and write them to their own raster for faster computing. This is helpful when comparing dates where a single Landsat date was excluded for analysis.

#landsat_test<-landsat[[n]]  ## The layer of interest in the landsat raster brick

#writeRaster(landsat_test, filename="./landsat_testday.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)

#starfm_test<- brick(starfm[[n]]) ## The layer of interest in the starfm raster brick

#writeRaster(starfm_test, filename="./starfm_testday.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)


### Set the landsat pixels with no data (value here set to -11111) to NA
#landsat<- landsat_test
landsat_masked<- Which(landsat == 0, cells = TRUE)
landsat[landsat_masked] <- -11111
landsat_plot<- reclassify(landsat,cbind(-Inf, -11111, NA))
landsat_plot<- landsat_plot/10000 ## Divide by 10000 to get NDSI data into the range of -1 to 1

### Here modis no data values are already set to NA. Only need to divide to get NDSI data into range.
modis_test<- modis[[8]]
modis_plot<- modis_test/10000

### Set the starfm pixels with no data (value here set to -11111) to NA
#starfm<- starfm_test
starfm_plot<- reclassify(starfm,cbind(-Inf, -11111, NA))
starfm_plot<- starfm_plot/10000

#### Plot the landsat NDSI values
#par(mfrow=c(1,2))
dev.new(height=0.91*nrow(landsat_plot)/50, width=1.09*ncol(landsat_plot)/50)
par(mar = c(5,5,5,3.5))
plot(landsat_plot, col = rev(cm.colors(40)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("NDSI from Landsat", cex = 2.5), line = 1.0)
#mtext("3 December 2015", line = 0.15, cex = 1.75)
plot(landsat_plot, legend.only = TRUE, col = rev(cm.colors(40)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", add = TRUE)


### Plot the STARFM NDSI values
#par(mar = c(3,3,4,3.5))
dev.new(height=0.91*nrow(starfm_plot)/50, width=1.09*ncol(starfm_plot)/50)
par(mar = c(5,5,5,3.5))
plot(starfm_plot, col = rev(cm.colors(40)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("NDSI from STARFM", cex = 2.5), line = 1.0)
#mtext("3 December 2015", line = 0.15, cex = 1.75)
plot(starfm_plot, legend.only = TRUE, col = rev(cm.colors(40)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", add = TRUE)

### Plot the modis NDSI values
par(mar = c(3,3,4,3.5))
plot(modis_plot, col = rev(cm.colors(40)), main = "NDSI from MODIS ")
mtext("Dec 3, 2015", line = 0.5)
plot(ER_proj, border = "black", add = TRUE)


## Find the fraction of data for each raster

## For landsat
area_L8<- landsat_plot@nrows*landsat_plot@ncols # total number of pixels for Landsat data
area_starfm<- starfm_plot@nrows*starfm_plot@ncols # total number of pixels for starfm data 

landsat_data<- landsat_plot
landsat_data[is.na(landsat_data)]<- 0  # set all NA values to zero
landsat_true<- Which(landsat_data != 0, cells = TRUE) 
landsat_data[landsat_true]<- 1  ## set all non-NA values to 1

landsat_valid<- as.vector(landsat_data, mode = 'numeric') # turn above data into vector for easy analysis 
L8_all_sum<- sum(landsat_valid) ## find the sum of all pixels with a value of 1

L8_data<- L8_all_sum/area_L8 ## find the fraction of pixels in the landsat raster that have data


## Same as above, but for STARFM raster
starfm_data<- starfm_plot
starfm_data[is.na(starfm_data)]<- 0
starfm_true<- Which(starfm_data != 0, cells = TRUE)
starfm_data[starfm_true]<- 1

starfm_valid<- as.vector(starfm_data, mode = 'numeric')
starfm_all_sum<- sum(starfm_valid)

model_data<- starfm_all_sum/area_starfm ## find the fraction of pixels in the starfm raster that have data


## Creating binary snow rasters

# Any values in our rasters of below 0.4 are given a new value of zero
landsat_binary<- reclassify(landsat_plot,cbind(-Inf, 0.4, 0))
starfm_binary<- reclassify(starfm_plot,cbind(-Inf, 0.4, 0))

# Any values in our rasters of 0.4 or above are given a new value of 1
landsat_binary<- reclassify(landsat_binary,cbind(0.4, Inf, 1))
starfm_binary<- reclassify(starfm_binary,cbind(0.4, Inf, 1))



## create levelplots to display the snow status. Requires the rasterVis library
landsat_bi_fac<- as.factor(landsat_binary)
rat<- levels(landsat_bi_fac)[[1]]
rat[["status"]] <- c("no snow", "snow")
levels(landsat_bi_fac)<- rat


dev.new()
levelplot(landsat_bi_fac, att = "status", col.regions = rev(cm.colors(40)), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), xlab = list(label = "Longitude", cex = 2), ylab = list(label = "Latitude", cex = 2),  main = list(label = "Landsat Snow Status", cex = 2),  colorkey = list(labels = list(height = 2, cex = 1.9))) + layer(sp.polygons(ER_proj)) 



starfm_ratify<- ratify(starfm_binary)
starfm_bi_fac<- as.factor(starfm_ratify)
rat<- levels(starfm_bi_fac)[[1]]
rat[["status"]] <- c("no snow", "snow")
levels(starfm_bi_fac)<- rat

dev.new()
levelplot(starfm_bi_fac, att = "status", col.regions = rev(cm.colors(40)), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), xlab = list(label = "Longitude", cex = 2), ylab = list(label = "Latitude", cex = 2), main = list(label = "STARFM Snow Status", cex = 2), colorkey = list(labels = list(height = 1, cex = 1.9))) + layer(sp.polygons(ER_proj))




## Find the fSCA for each raster

L8_snow<- as.vector(landsat_binary, mode = 'numeric')
L8_snow_sum<- sum(L8_snow, na.rm = TRUE) ## Find the sum all pixels that have a value equal to 1 (snow status)

L8_fsca<- L8_snow_sum/area_L8  #Landsat fSCA


starfm_snow<- as.vector(starfm_binary, mode = 'numeric') ## Find the sum all pixels that have a value equal to 1 (snow status)
starfm_snow_sum<- sum(starfm_snow, na.rm = TRUE)

starfm_fsca<- starfm_snow_sum/area_starfm #STARFM fSCA



## Create a confusion matrix to evaluate the performance of the model against the Landsat data

landsat_factor<- as.factor(L8_snow)
starfm_factor<- as.factor(starfm_snow)

results<- confusionMatrix(starfm_factor,landsat_factor, positive = "1")

# Save the results in matrices that can easily be viewed for reference
overall<- as.matrix(results, what = "overall")
classes<- as.matrix(results, what = "classes")


# Want to know the individual values, so save these values as their own variable
acc_model<- overall[[1]]
spec_model<- classes[[2]]
prec_model<- classes[[5]]
recall_model<- classes[[6]]
f1_model<- classes[[7]]

date<- "12.3.2015"

results_end<- c(date, acc_model, spec_model, recall_model, prec_model, f1_model, L8_fsca, starfm_fsca, L8_data, model_data)
labels<- c("Date", "Accuracy", "Specificity", "Recall", "Precision", "F-Score", "L8 fSCA", "STARFM fSCA", "L8 Data", "STARFM Data")

all_data<- data.frame(labels, results_end)
