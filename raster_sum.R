setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/model_id')


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


############# Put all the rasters of the id'd dates (rasters of 1 for correct prediction, 0 for incorrect prediction), into one brick raster

#brick<- do.call(brick, lapply(list.files(path = getwd(), pattern = "id_*"), raster))

## Sum all the rasters together, so the value of the sum for each pixel is the value of the number of instances that the model prediction was correct

#all_sum<- calc(brick, fun = sum, na.rm = TRUE)

## Load in the other data that will be needed for analysis (watershed shapefile, DEM, calculated aspect, calculated slope)

east_DEM<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/East_DEM.tif')

east_slope<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/east_slope.tif')

east_aspect<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/east_aspect.tif')
  
ER<- readOGR('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/ER_watershed_shp/EastRiver_Project.shp') #shapefile of the study watershed for reference

## Sum of correct instances
east_sum<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/east_sum.tif')

## Sum of total number of predictions, regardless of correctness
east_model_sum<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/model_id/model_pred_sum.tif')


## Project the shapefile to the same projection as the rasters

data_proj<- crs(east_sum)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

## Plot the value of the pixels sums and write it to it's own raster to use later

dev.new(height=0.91*nrow(east_sum)/50, width=1.09*ncol(east_sum)/50)
par(mar = c(5,5,5,3.5))
plot(east_sum, col = brewer.pal(9, name = "OrRd"), main = "", xlab = "Longitude", ylab = "Longitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Sums of Correct Pixels", cex = 2.5), line = 1.0)
plot(east_sum, legend.only = TRUE, col = brewer.pal(9, name = "OrRd"), axis.arg = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


dev.new(height=0.91*nrow(east_model_sum)/50, width=1.09*ncol(east_model_sum)/50)
par(mar = c(5,5,5,3.5))
plot(east_model_sum, col = brewer.pal(9, name = "OrRd"), main = "", xlab = "Longitude", ylab = "Longitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Total Predictions per Pixel", cex = 2.5), line = 1.0)
plot(east_model_sum, legend.only = TRUE, col = brewer.pal(9, name = "OrRd"), axis.arg = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


#writeRaster(all_sum, filename="C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/DEM_Analysis/east_sum.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)

########### Find the difference between the number of times a model prediction was able to be made on a pixel and the number of times that model prediction was correct

sum_diff<- east_model_sum - east_sum

dev.new(height=0.91*nrow(sum_diff)/50, width=1.09*ncol(sum_diff)/50)
par(mar = c(5,5,5,3.5))
plot(sum_diff, col = brewer.pal(9, name = "OrRd"),  main = "", xlab = "Longitude", ylab = "Longitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Difference Total and Correct Predictions", cex = 1.5), line = 1.0)
plot(sum_diff, legend.only = TRUE, col = brewer.pal(9, name = "OrRd"), axis.arg = list(cex.axis = 1.9))
plot(ER_proj, border = 'black', lwd = 2, add = TRUE)


perc_correct<- east_sum/east_model_sum

dev.new(height=0.91*nrow(perc_correct)/50, width=1.09*ncol(perc_correct)/50)
par(mar = c(5,5,5,3.5))
plot(perc_correct, main = "", xlab = "Longitude", ylab = "Longitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy by Pixel", cex = 2), line = 1.0)
plot(perc_correct, legend.only = TRUE, axis.arg = list(cex.axis = 1.9))
plot(ER_proj, border = 'black', lwd = 2, add = TRUE)



########### Analyses by elevation

plot(east_DEM, main = "Elevation (m)")
plot(ER_proj, border = 'black', add = TRUE)

## Create elevation bands by dividing the data into thirds

# Find the min and max elevations

 elev_stats<- stats(east_DEM)
 elev_min<- elev_stats[[4]]
 elev_max<- elev_stats[[8]]
 
 elev_range<- elev_max-elev_min
# elev_cat<- elev_range/3
# low_value <- elev_min + elev_cat
# med_value <- low_value + elev_cat

## Create a hypsometric curve for elevations
east_spg<- as(east_DEM, "SpatialGridDataFrame")
hypsometric(east_spg)


## Calculate the empirical cumulative distribution of the elevation data

DEM_vec<- as.vector(east_DEM)
ecdf_DEM<- ecdf(DEM_vec)
plot(ecdf_DEM, main = "ECDF of East River Elevations", xlab = "Meters")

y<- ecdf(DEM_vec);quantile(y, c(0.33, 0.66))



## To create the ecdf another way

x<- sort(DEM_vec)

e_cdf<- 1:length(x)/ length(x)
#plot(x, e_cdf, type = 's')

low_elev_limit<- x[which(e_cdf >= 0.33)[1]]
mid_elev_limit<- x[which(e_cdf >= 0.66)[1]]



# Create the elevation bands with the ecdf values

elev_low<- (east_DEM <= low_elev_limit)
elev_med<- (east_DEM > low_elev_limit & east_DEM <= mid_elev_limit)
elev_high<- (east_DEM > mid_elev_limit)

elev_bands <- east_DEM
elev_bands[elev_low] <- 1
elev_bands[elev_med] <- 2
elev_bands[elev_high] <- 3

plot(elev_bands, main = "Elevation Bands")
plot(ER_proj, border = 'black', add = TRUE)

elev_bands_fac<- as.factor(elev_bands)
rat<- levels(elev_bands_fac)[[1]] 
rat[["bands"]] <- c("low", "medium", "high")
levels(elev_bands_fac)<- rat

dev.new()
levelplot(elev_bands_fac, col.regions = rev(terrain.colors(3)), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), xlab = list(label = "Longitude", cex = 2), ylab = list(label = "Latitude", cex = 2), main = list(label = "Elevation Bands", cex = 2), colorkey = list(labels = list(height = 1, cex = 1.7))) + layer(sp.polygons(ER_proj, lwd = 2))


## Plot the number of correct model instances by elevation band

# For low elevations

low_sums <- east_sum
masked <- Which(elev_low == 0, cells=TRUE)
low_sums[masked] <- NA

plot(low_sums, main = "Sum values at Low Elev")
plot(ER_proj, border = 'black', add = TRUE)


elev_low_perc <- perc_correct
masked <- Which(elev_low == 0, cells=TRUE)
elev_low_perc[masked] <- NA


dev.new(height=0.91*nrow(elev_low_perc)/50, width=1.09*ncol(elev_low_perc)/50)
par(mar = c(5,5,5,3.5))
plot(elev_low_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at Low Elev", cex = 2.5), line = 1.0)
plot(elev_low_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For mid elevations

# med_sums <- east_sum
# masked <- Which(elev_med == 0, cells=TRUE)
# med_sums[masked] <- NA
# 
# plot(med_sums, main = "Sum values at Med Elev")
# plot(ER_proj, border = 'black', add = TRUE)


elev_med_perc <- perc_correct
masked <- Which(elev_med == 0, cells=TRUE)
elev_med_perc[masked] <- NA


dev.new(height=0.91*nrow(elev_med_perc)/50, width=1.09*ncol(elev_med_perc)/50)
par(mar = c(5,5,5,3.5))
plot(elev_med_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at Mid Elev", cex = 2.5), line = 1.0)
plot(elev_med_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For high elevations

# high_sums <- east_sum
# masked <- Which(elev_high == 0, cells=TRUE)
# high_sums[masked] <- NA
# 
# plot(high_sums, main = "Sum values at High Elev")
# plot(ER_proj, border = 'black', add = TRUE)


elev_high_perc <- perc_correct
masked <- Which(elev_high == 0, cells=TRUE)
elev_high_perc[masked] <- NA


dev.new(height=0.91*nrow(elev_med_perc)/50, width=1.09*ncol(elev_med_perc)/50)
par(mar = c(5,5,5,3.5))
plot(elev_high_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at High Elev", cex = 2.5), line = 1.0)
plot(elev_high_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### Create box and whisker plots to better analyze data availability by elevation band

low_elev_data <- east_model_sum
masked <- Which(elev_low == 0, cells=TRUE)
low_elev_data[masked] <- NA


dev.new(height=0.91*nrow(low_elev_data)/50, width=1.09*ncol(low_elev_data)/50)
par(mar = c(5,5,5,3.5))
plot(low_elev_data, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at Low Elev", cex = 2.5), line = 1.0)
plot(low_elev_data, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



low_elev_data_vec<- as.vector(low_elev_data)

boxplot(low_elev_data_vec, main = "Data available at low elevation", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE, notch = TRUE)


mid_elev_data <- east_model_sum
masked <- Which(elev_med == 0, cells=TRUE)
mid_elev_data[masked] <- NA


dev.new(height=0.91*nrow(mid_elev_data)/50, width=1.09*ncol(mid_elev_data)/50)
par(mar = c(5,5,5,3.5))
plot(mid_elev_data,  col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at Mid Elev", cex = 2.5), line = 1.0)
plot(mid_elev_data, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


mid_elev_data_vec<- as.vector(mid_elev_data)

boxplot(mid_elev_data_vec, main = "Data available at mid elevation", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE, notch = TRUE)


high_elev_data <- east_model_sum
masked <- Which(elev_high == 0, cells=TRUE)
high_elev_data[masked] <- NA


dev.new(height=0.91*nrow(high_elev_data)/50, width=1.09*ncol(high_elev_data)/50)
par(mar = c(5,5,5,3.5))
plot(high_elev_data, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at High Elev", cex = 2.5), line = 1.0)
plot(high_elev_data, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


high_elev_data_vec<- as.vector(high_elev_data)

boxplot(high_elev_data_vec, main = "Data available at high elevation", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE)


boxplot(low_elev_data_vec, mid_elev_data_vec, high_elev_data_vec, 
        main = "Data Available at Elevation Bands", cex.main = 1.8,
        at = c(1,2,3),
        names = c('Low', 'Mid', 'High'),
        las = 2,
        xlab = "# of Predictions", cex.lab = 1.5, cex.axis = 1.5, 
        col = terrain.colors(3), 
        horizontal = TRUE)




############# Analyses by slope

 slope_stats<- stats(east_slope)
 slope_max<- slope_stats[[8]]
 slope_min<- slope_stats[[4]]
 
 slope_cat<- slope_max/3
# med_deg <- slope_min + slope_cat*2

## Calculate the empirical cumulative distribution of the elevation data

slope_vec<- as.vector(east_slope)
ecdf_slope<- ecdf(slope_vec)
plot(ecdf_slope, main = "ECDF of East River Slopes", xlab = "Degrees")

y<- ecdf(slope_vec);quantile(y, c(0.33, 0.66))



## To create the ecdf another way

x<- sort(slope_vec)

e_cdf<- 1:length(x)/ length(x)
#plot(x, e_cdf, type = 's')

low_slope_limit<- x[which(e_cdf >= 0.33)[1]]
mid_slope_limit<- x[which(e_cdf >= 0.66)[1]]





# Create the slope bands
slope_low<- (east_slope <= low_slope_limit)
slope_med<- (east_slope > low_slope_limit & east_slope <= mid_slope_limit)
slope_high<- (east_slope > mid_slope_limit)

slope_bands <- east_slope
slope_bands[slope_low] <- 1
slope_bands[slope_med] <- 2
slope_bands[slope_high] <- 3

plot(slope_bands, main = "Slope Bands")
plot(ER_proj, border = 'black', add = TRUE)

slope_bands_fac<- as.factor(slope_bands)
rat<- levels(slope_bands_fac)[[1]] 
rat[["bands"]] <- c("low", "medium", "high")
levels(slope_bands_fac)<- rat

dev.new()
levelplot(slope_bands_fac, col.regions = rev(brewer.pal(3, name = "RdYlBu")), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), xlab = list(label = "Longitude", cex = 2), ylab = list(label = "Latitude", cex = 2), main = list(label = "Slope Bands", cex = 2), colorkey = list(labels = list(height = 1, cex = 1.7))) + layer(sp.polygons(ER_proj, lwd = 2))


## Plot the number of correct model instances by slope band

# For low slopes

# low_sums <- east_sum
# masked <- Which(slope_low == 0, cells=TRUE)
# low_sums[masked] <- NA
# 
# low_sums_crop<- crop(low_sums, extent(low_sums,2, 1590, 2, 1805))
# 
# plot(low_sums_crop, main = "Sum values at Low Slopes")
# plot(ER_proj, border = 'black', add = TRUE)


slope_low_perc <- perc_correct
masked <- Which(slope_low == 0, cells=TRUE)
slope_low_perc[masked] <- NA

slope_low_perc_crop<- crop(slope_low_perc, extent(slope_low_perc,2, 1590, 2, 1805))


dev.new(height=0.91*nrow(slope_low_perc_crop)/50, width=1.09*ncol(slope_low_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(slope_low_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at Low Slope", cex = 2.5), line = 1.0)
plot(slope_low_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



# For mid slopes

# med_sums <- east_sum
# masked <- Which(slope_med == 0, cells=TRUE)
# med_sums[masked] <- NA
# 
# med_sums_crop<- crop(med_sums, extent(med_sums,2, 1590, 2, 1805))
# 
# plot(med_sums_crop, main = "Sum values at Med Slopes")
# plot(ER_proj, border = 'black', add = TRUE)


slope_mid_perc <- perc_correct
masked <- Which(slope_med == 0, cells=TRUE)
slope_mid_perc[masked] <- NA

slope_mid_perc_crop<- crop(slope_mid_perc, extent(slope_mid_perc,2, 1590, 2, 1805))


dev.new(height=0.91*nrow(slope_mid_perc_crop)/50, width=1.09*ncol(slope_mid_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(slope_mid_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at Mid Slope", cex = 2.5), line = 1.0)
plot(slope_mid_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



# For high slopes

# high_sums <- east_sum
# masked <- Which(slope_high == 0, cells=TRUE)
# high_sums[masked] <- NA
# 
# high_sums_crop<- crop(high_sums, extent(high_sums,2, 1590, 2, 1805))
# 
# plot(high_sums_crop, main = "Sum values at High Slope")
# plot(ER_proj, border = 'black', add = TRUE)


slope_high_perc <- perc_correct
masked <- Which(slope_high == 0, cells=TRUE)
slope_high_perc[masked] <- NA

slope_high_perc_crop<- crop(slope_high_perc, extent(slope_high_perc,2, 1590, 2, 1805))


dev.new(height=0.91*nrow(slope_high_perc_crop)/50, width=1.09*ncol(slope_high_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(slope_high_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at High Slope", cex = 2.5), line = 1.0)
plot(slope_high_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



########### Analyses by aspect

plot(east_aspect, main = "Aspect (deg)")
plot(ER_proj, border = 'black', add = TRUE)

## Assign directions to variables

aspect_stats<- stats(east_aspect)

# north<- (east_aspect >= 0 & east_aspect <= 22.5 | east_aspect > 337.5)
# northeast<- (east_aspect > 22.5 & east_aspect <= 67.5)
# east<- (east_aspect > 67.5 & east_aspect <= 112.5)
# southeast<- (east_aspect > 112.5 & east_aspect <= 157.5)
# south<- (east_aspect > 157.5 & east_aspect <= 202.5)
# southwest<- (east_aspect > 202.5 & east_aspect <= 247.5)
# west<- (east_aspect > 247.5 & east_aspect <= 292.5)
# northwest<- (east_aspect > 292.5 & east_aspect <= 337.5)

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

dev.new()
levelplot(aspect_cat, att = 'direction', col.regions = c("#FF0000", "#FFA500", "#F0E68C", "#87CEEB"), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), xlab = list(label = "Longitude", cex = 2), ylab = list(label = "Latitude", cex = 2), main = list(label = 'Aspect', cex = 2), colorkey = list(labels = list(height = 1, cex = 1.7))) + layer(sp.polygons(ER_proj, lwd = 2))



## Plot the number of correct model instances by direction

# For north facing slopes

# north_sums <- east_sum
# masked <- Which(north == 0, cells=TRUE)
# north_sums[masked] <- NA
# 
# north_sums_crop<- crop(north_sums, extent(north_sums,2, 1590, 2, 1805))
# 
# plot(north_sums_crop, main = "Sum values at North facing Slopes")
# plot(ER_proj, border = 'black', add = TRUE)


aspect_north_perc <- perc_correct
masked <- Which(north == 0, cells=TRUE)
aspect_north_perc[masked] <- NA

aspect_north_perc_crop<- crop(aspect_north_perc, extent(aspect_north_perc,2, 1590, 2, 1805))


dev.new(height=0.91*nrow(aspect_north_perc_crop)/50, width=1.09*ncol(aspect_north_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(aspect_north_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at North Aspect", cex = 2.5), line = 1.0)
plot(aspect_north_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



# For east facing slopes

# east_sums <- east_sum
# masked <- Which(east == 0, cells=TRUE)
# east_sums[masked] <- NA
# 
# east_sums_crop<- crop(east_sums, extent(east_sums,2, 1590, 2, 1805))
# 
# plot(east_sums_crop, main = "Sum values at East facing Slope")
# plot(ER_proj, border = 'black', add = TRUE)


aspect_east_perc <- perc_correct
masked <- Which(east == 0, cells=TRUE)
aspect_east_perc[masked] <- NA

aspect_east_perc_crop<- crop(aspect_east_perc, extent(aspect_east_perc,2, 1590, 2, 1805))

dev.new(height=0.91*nrow(aspect_east_perc_crop)/50, width=1.09*ncol(aspect_east_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(aspect_east_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at East Aspect", cex = 2.5), line = 1.0)
plot(aspect_east_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For south facing slopes

# south_sums <- east_sum
# masked <- Which(south == 0, cells=TRUE)
# south_sums[masked] <- NA
# 
# south_sums_crop<- crop(south_sums, extent(south_sums,2, 1590, 2, 1805))
# 
# plot(south_sums_crop, main = "Sum values at South facing Slope")
# plot(ER_proj, border = 'black', add = TRUE)

aspect_south_perc <- perc_correct
masked <- Which(south == 0, cells=TRUE)
aspect_south_perc[masked] <- NA

aspect_south_perc_crop<- crop(aspect_south_perc, extent(aspect_south_perc,2, 1590, 2, 1805))

dev.new(height=0.91*nrow(aspect_south_perc_crop)/50, width=1.09*ncol(aspect_south_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(aspect_south_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at South Aspect", cex = 2.5), line = 1.0)
plot(aspect_south_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For west facing slopes

# west_sums <- east_sum
# masked <- Which(west == 0, cells=TRUE)
# west_sums[masked] <- NA
# 
# west_sums_crop<- crop(west_sums, extent(west_sums,2, 1590, 2, 1805))
# 
# plot(west_sums_crop, main = "Sum values at West facing Slope")
# plot(ER_proj, border = 'black', add = TRUE)


aspect_west_perc <- perc_correct
masked <- Which(west == 0, cells=TRUE)
aspect_west_perc[masked] <- NA

aspect_west_perc_crop<- crop(aspect_west_perc, extent(aspect_west_perc,2, 1590, 2, 1805))

dev.new(height=0.91*nrow(aspect_west_perc_crop)/50, width=1.09*ncol(aspect_west_perc_crop)/50)
par(mar = c(5,5,5,3.5))
plot(aspect_west_perc_crop, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy at West Aspect", cex = 2.5), line = 1.0)
plot(aspect_west_perc_crop, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### Create box and whisker plots to better analyze data availability by aspect direction

north_data <- east_model_sum
masked <- Which(north == 0, cells=TRUE)
north_data[masked] <- NA

north_data_crop<- crop(north_data, extent(north_data, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(north_data_crop)/50, width=1.09*ncol(north_data_crop)/50)
par(mar = c(5,5,5,3.5))
plot(north_data_crop, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at North Aspect", cex = 2.5), line = 1.0)
plot(north_data_crop, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


north_data_vec<- as.vector(north_data_crop)

boxplot(north_data_vec, main = "Data available at north aspect", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE, notch = TRUE)


east_data <- east_model_sum
masked <- Which(east == 0, cells=TRUE)
east_data[masked] <- NA

east_data_crop<- crop(east_data, extent(east_data, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(east_data_crop)/50, width=1.09*ncol(east_data_crop)/50)
par(mar = c(5,5,5,3.5))
plot(east_data_crop, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at East Aspect", cex = 2.5), line = 1.0)
plot(east_data_crop, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


east_data_vec<- as.vector(east_data_crop)

boxplot(east_data_vec, main = "Data available at east aspect", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE, notch = TRUE)


south_data <- east_model_sum
masked <- Which(south == 0, cells=TRUE)
south_data[masked] <- NA

south_data_crop<- crop(south_data, extent(south_data, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(south_data_crop)/50, width=1.09*ncol(south_data_crop)/50)
par(mar = c(5,5,5,3.5))
plot(south_data_crop, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at South Aspect", cex = 2.5), line = 1.0)
plot(south_data_crop, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


south_data_vec<- as.vector(south_data_crop)

boxplot(south_data_vec, main = "Data available at south aspect", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE)


west_data <- east_model_sum
masked <- Which(west == 0, cells=TRUE)
west_data[masked] <- NA

west_data_crop<- crop(west_data, extent(west_data, 2, 1590, 2, 1805))

dev.new(height=0.91*nrow(west_data_crop)/50, width=1.09*ncol(west_data_crop)/50)
par(mar = c(5,5,5,3.5))
plot(west_data_crop, col = rev(plasma(25)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Data Available at West Aspect", cex = 2.5), line = 1.0)
plot(west_data_crop, legend.only = TRUE, col = rev(plasma(25)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


west_data_vec<- as.vector(west_data_crop)

boxplot(west_data_vec, main = "Data available at west aspect", xlab = "# of Predictions", col = 'orange', border = 'brown', horizontal = TRUE)


boxplot(north_data_vec, east_data_vec, south_data_vec, west_data_vec, 
        main = "Data Available at Aspect", cex.main = 1.8,
        at = c(1,2,3,4),
        names = c('North', 'East', 'South', 'West'),
        las = 2,
        xlab = "# of Predictions", cex.lab = 1.5, cex.axis = 1.4,
        col = c("#FF0000", "#FFA500", "#F0E68C", "#87CEEB"), 
        horizontal = TRUE)



############ Create barplots to show the means of correct predictions by spatial feature

## Convert all values from rasters into vectors for analysis

low_elev_vec<- as.vector(elev_low_perc)
mid_elev_vec<- as.vector(elev_med_perc)
high_elev_vec<- as.vector(elev_high_perc)


low_slope_vec<- as.vector(slope_low_perc_crop)
mid_slope_vec<- as.vector(slope_mid_perc_crop)
high_slope_vec<- as.vector(slope_high_perc_crop)

north_aspect_vec<- as.vector(aspect_north_perc_crop)
east_aspect_vec<- as.vector(aspect_east_perc_crop)
south_aspect_vec<- as.vector(aspect_south_perc_crop)
west_aspect_vec<- as.vector(aspect_west_perc_crop)


## Find the means and standard devs for each landscape characteristic

low_elev_mean<- mean(low_elev_vec, na.rm = TRUE)
low_elev_sd<- sd(low_elev_vec, na.rm = TRUE)
mid_elev_mean<- mean(mid_elev_vec, na.rm = TRUE)
mid_elev_sd<- sd(mid_elev_vec, na.rm = TRUE)
high_elev_mean<- mean(high_elev_vec, na.rm = TRUE)
high_elev_sd<- sd(high_elev_vec, na.rm = TRUE)

low_slope_mean<- mean(low_slope_vec, na.rm = TRUE)
low_slope_sd<- sd(low_slope_vec, na.rm = TRUE)
mid_slope_mean<- mean(mid_slope_vec, na.rm = TRUE)
mid_slope_sd<- sd(mid_slope_vec, na.rm = TRUE)
high_slope_mean<- mean(high_slope_vec, na.rm = TRUE)
high_slope_sd<- sd(high_slope_vec, na.rm = TRUE)

north_aspect_mean<- mean(north_aspect_vec, na.rm = TRUE)
north_aspect_sd<- sd(north_aspect_vec, na.rm = TRUE)
east_aspect_mean<- mean(east_aspect_vec, na.rm = TRUE)
east_aspect_sd<- sd(east_aspect_vec, na.rm = TRUE)
south_aspect_mean<- mean(south_aspect_vec, na.rm = TRUE)
south_aspect_sd<- sd(south_aspect_vec, na.rm = TRUE)
west_aspect_mean<- mean(west_aspect_vec, na.rm = TRUE)
west_aspect_sd<- sd(west_aspect_vec, na.rm = TRUE)

############ Create barplots to show the means of correct predictions by spatial feature

## Convert all values from rasters into vectors for analysis

low_elev_vec<- as.vector(elev_low_perc)
mid_elev_vec<- as.vector(elev_med_perc)
high_elev_vec<- as.vector(elev_high_perc)


low_slope_vec<- as.vector(slope_low_perc_crop)
mid_slope_vec<- as.vector(slope_mid_perc_crop)
high_slope_vec<- as.vector(slope_high_perc_crop)

north_aspect_vec<- as.vector(aspect_north_perc_crop)
east_aspect_vec<- as.vector(aspect_east_perc_crop)
south_aspect_vec<- as.vector(aspect_south_perc_crop)
west_aspect_vec<- as.vector(aspect_west_perc_crop)


## Find the means and standard devs for each landscape characteristic

low_elev_mean<- mean(low_elev_vec, na.rm = TRUE)
low_elev_sd<- sd(low_elev_vec, na.rm = TRUE)
mid_elev_mean<- mean(mid_elev_vec, na.rm = TRUE)
mid_elev_sd<- sd(mid_elev_vec, na.rm = TRUE)
high_elev_mean<- mean(high_elev_vec, na.rm = TRUE)
high_elev_sd<- sd(high_elev_vec, na.rm = TRUE)

low_slope_mean<- mean(low_slope_vec, na.rm = TRUE)
low_slope_sd<- sd(low_slope_vec, na.rm = TRUE)
mid_slope_mean<- mean(mid_slope_vec, na.rm = TRUE)
mid_slope_sd<- sd(mid_slope_vec, na.rm = TRUE)
high_slope_mean<- mean(high_slope_vec, na.rm = TRUE)
high_slope_sd<- sd(high_slope_vec, na.rm = TRUE)

north_aspect_mean<- mean(north_aspect_vec, na.rm = TRUE)
north_aspect_sd<- sd(north_aspect_vec, na.rm = TRUE)
east_aspect_mean<- mean(east_aspect_vec, na.rm = TRUE)
east_aspect_sd<- sd(east_aspect_vec, na.rm = TRUE)
south_aspect_mean<- mean(south_aspect_vec, na.rm = TRUE)
south_aspect_sd<- sd(south_aspect_vec, na.rm = TRUE)
west_aspect_mean<- mean(south_aspect_vec, na.rm = TRUE)
west_aspect_sd<- sd(south_aspect_vec, na.rm = TRUE)

## Create and display the bar plots

## Elevation bar plot


labels<- c("Low", 'Mid', 'High')
means<- c(low_elev_mean, mid_elev_mean, high_elev_mean)
sds<- c(low_elev_sd, mid_elev_sd, high_elev_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg = labels, cex.names = 1.2, main="Accuracies by Elevation Band", col = terrain.colors(3), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.5, las=1, ylim=range(pretty(c(0, 1.1))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 0.15, paste(round(means, digits = 2)), cex = 1.2)


## Slope bar plot

labels<- c("Low", 'Mid', 'High')
means<- c(low_slope_mean, mid_slope_mean, high_slope_mean)
sds<- c(low_slope_sd, mid_slope_sd, high_slope_sd)


par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Accuracies by Slope Band", col = rev(brewer.pal(3, name = "RdYlBu")), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.5, las=1, ylim=range(pretty(c(0, 1.1))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 0.2, paste(round(means, digits = 2)), cex = 1.2)

## Aspect bar plot

labels<- c("North", 'East', 'South', 'West')
means<- c(north_aspect_mean, east_aspect_mean, south_aspect_mean, west_aspect_mean)
sds<- c(north_aspect_sd, east_aspect_sd, south_aspect_sd, west_aspect_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Accuracies by Aspect", col = c("#FF0000", "#FFA500", "#F0E68C", "#87CEEB"), ylab="Mean Value", cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.5, las=1, ylim= range(pretty(c(0, 1.1))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 0.2, paste(round(means, digits = 2)), cex = 1.2)



########### Analysis by vegetation class



east_veg_class<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/US_200EVT_ERarea/LF_fullarea/Lf_fullarea.tif')

#east_veg_class<- projectRaster(east_veg_class, crs = east_sum)

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



### Now we can map the percent correct by vegetation class


east_trees <- (east_veg_class == 1)
east_shrub <- (east_veg_class == 2)
east_grass <- (east_veg_class == 3 | east_veg_class == 8)
east_clear <- (east_veg_class == 5 | east_veg_class == 9)
east_snow_ice <- (east_veg_class == 7)
east_other<- (east_veg_class == 4 | east_veg_class == 6 | east_veg_class == 7)


## Same as above, just all in one raster
veg_class<- east_veg_class
veg_class[east_trees] <- 1
veg_class[east_shrub] <- 2
veg_class[east_grass] <- 3
veg_class[east_clear] <- 4
veg_class[east_other]<- 5



veg_class_fac<- as.factor(veg_class)
rat_veg<- levels(veg_class_fac)[[1]] 
rat_veg[["class"]] <- c("trees", "shrubs", "grass", "clear", "other")
levels(veg_class_fac)<- rat_veg

dev.new()
levelplot(veg_class_fac,  col.regions = c("#33A02C", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#6A3D9A"), scales = list(x = list(cex = 1.5), y = list(cex = 1.5, rot = 90)), xlab = list(label = "Meters", cex = 2), ylab = list(label = "Meters", cex = 2), main = list(label = "Landcover Classes", cex = 2), colorkey = list(labels = list(height = 1, cex = 1.5))) + layer(sp.polygons(ER_veg, lwd = 2))




trees_proj<- projectRaster(east_trees, crs = data_proj, method = 'ngb')
trees_proj <- resample(trees_proj, perc_correct, method = 'ngb')

shrubs_proj<- projectRaster(east_shrub, crs = data_proj, method = 'ngb')
shrubs_proj <- resample(shrubs_proj, perc_correct, method = 'ngb')

grass_proj<- projectRaster(east_grass, crs = data_proj, method = 'ngb')
grass_proj <- resample(grass_proj, perc_correct, method = 'ngb')

clear_proj<- projectRaster(east_clear, crs = data_proj, method = 'ngb')
clear_proj <- resample(clear_proj, perc_correct, method = 'ngb')

snow_ice_proj<- projectRaster(east_snow_ice, crs = data_proj, method = 'ngb')
snow_ice_proj <- resample(snow_ice_proj, perc_correct, method = 'ngb')

other_proj<- projectRaster(east_other, crs = data_proj, method = 'ngb')
other_proj <- resample(other_proj, perc_correct, method = 'ngb')


trees_perc <- perc_correct
masked <- Which(trees_proj == 0, cells=TRUE)
trees_perc[masked] <- NA

dev.new(height=0.91*nrow(trees_perc)/50, width=1.09*ncol(trees_perc)/50)
par(mar = c(5,5,5,3.5))
plot(trees_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy for Tree Areas", cex = 2.5), line = 1.0)
plot(trees_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


shrubs_perc <- perc_correct
masked <- Which(shrubs_proj == 0, cells=TRUE)
shrubs_perc[masked] <- NA


dev.new(height=0.91*nrow(shrubs_perc)/50, width=1.09*ncol(shrubs_perc)/50)
par(mar = c(5,5,5,3.5))
plot(shrubs_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy for Shrub Areas", cex = 2.5), line = 1.0)
plot(shrubs_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


grass_perc <- perc_correct
masked <- Which(grass_proj == 0, cells=TRUE)
grass_perc[masked] <- NA


dev.new(height=0.91*nrow(grass_perc)/50, width=1.09*ncol(grass_perc)/50)
par(mar = c(5,5,5,3.5))
plot(grass_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy for Grass Areas", cex = 2.5), line = 1.0)
plot(grass_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


clear_perc <- perc_correct
masked <- Which(clear_proj == 0, cells=TRUE)
clear_perc[masked] <- NA


dev.new(height=0.91*nrow(clear_perc)/50, width=1.09*ncol(clear_perc)/50)
par(mar = c(5,5,5,3.5))
plot(clear_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy for Clear Areas", cex = 2.5), line = 1.0)
plot(clear_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


snow_ice_perc <- perc_correct
masked <- Which(snow_ice_proj == 0, cells=TRUE)
snow_ice_perc[masked] <- NA
 
plot(snow_ice_perc, main = "Perc correct for Snow/Ice Areas")
plot(ER_proj, border = 'black', add = TRUE)


other_perc <- perc_correct
masked <- Which(other_proj == 0, cells=TRUE)
other_perc[masked] <- NA


dev.new(height=0.91*nrow(other_perc)/50, width=1.09*ncol(other_perc)/50)
par(mar = c(5,5,5,3.5))
plot(other_perc, main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("Accuracy for All Other Areas", cex = 2.5), line = 1.0)
plot(other_perc, legend.only = TRUE, axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

############ Create barplots to show the means of correct predictions by spatial feature

## Convert all values from rasters into vectors for analysis

trees_perc_vec<- as.vector(trees_perc)
shrubs_perc_vec<- as.vector(shrubs_perc)
grass_perc_vec<- as.vector(grass_perc)
clear_perc_vec<- as.vector(clear_perc)
snow_ice_perc_vec<- as.vector(snow_ice_perc)
other_perc_vec<- as.vector(other_perc)


## Find the means and standard devs for each landscape characteristic

trees_mean<- mean(trees_perc_vec, na.rm = TRUE)
trees_sd<- sd(trees_perc_vec, na.rm = TRUE)
shrubs_mean<- mean(shrubs_perc_vec, na.rm = TRUE)
shrubs_sd<- sd(shrubs_perc_vec, na.rm = TRUE)
grass_mean<- mean(grass_perc_vec, na.rm = TRUE)
grass_sd<- sd(grass_perc_vec, na.rm = TRUE)
clear_mean<- mean(clear_perc_vec, na.rm = TRUE)
clear_sd<- sd(clear_perc_vec, na.rm = TRUE)
snow_ice_mean<- mean(snow_ice_perc_vec, na.rm = TRUE)
snow_ice_sd<- sd(snow_ice_perc_vec, na.rm = TRUE)
other_mean<- mean(other_perc_vec, na.rm = TRUE)
other_sd<- sd(other_perc_vec, na.rm = TRUE)


## Vegetation accuracies bar plot

labels<- c('Trees', 'Shrubs', 'Grass', 'Clear', 'Other')
means<- c(trees_mean, shrubs_mean, grass_mean, clear_mean, other_mean)
sds<- c(trees_sd, shrubs_sd, grass_sd, clear_sd, other_sd)

par(mfrow=c(1,1))
mids<- barplot(means, names.arg=labels, cex.names = 1.2, main="Accuracies by Landcover Class", ylab="Mean Value", col = c("#33A02C", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#6A3D9A"), cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.5, las=1, ylim= range(pretty(c(0, 1.1))))
arrows(x0=mids, y0=means-sds, x1=mids, y1=means+sds, code=3, angle=90, length=0.1)
text(mids, means + 0.15, paste(round(means, digits = 2)), cex = 1.2)


## Find the percentage of area that each vegetation type covers

total_area<- ncol(east_sum) * nrow(east_sum)

trees_vec<- as.vector(trees_proj)
shrubs_vec<- as.vector(shrubs_proj)
grass_vec<- as.vector(grass_proj)
clear_vec<- as.vector(clear_proj)
snow_ice_vec<- as.vector(snow_ice_proj)
other_vec<- as.vector(other_proj)

trees_sum<- sum(trees_vec, na.rm = TRUE)
shrubs_sum<- sum(shrubs_vec, na.rm = TRUE)
grass_sum<- sum(grass_vec, na.rm = TRUE)
clear_sum<- sum(clear_vec, na.rm = TRUE)
snow_ice_sum<- sum(snow_ice_vec, na.rm = TRUE)
other_sum<- sum(other_vec, na.rm = TRUE)

trees_area<- (trees_sum/total_area) *100
shrubs_area<- (shrubs_sum/total_area) *100
grass_area<- (grass_sum/total_area) *100
clear_area<- (clear_sum/total_area) * 100
snow_ice_area<- snow_ice_sum/total_area
other_area<- (other_sum/total_area) *100


labels<- c('Trees', 'Shrubs', 'Grass', 'Clear', 'Other')
areas<- c(trees_area, shrubs_area, grass_area, clear_area, other_area)

par(mfrow=c(1,1))
mids<- barplot(areas, names.arg= labels, main="Percent Area by Landcover Class", ylab="Percent", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.1, col = c("#33A02C", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#6A3D9A"), las=1, ylim= range(pretty(c(0, 100))))
text(mids, areas + 7, paste(round(areas, digits = 0)), cex = 1.5)


### Find and plot data availability by landcover type

tree_data <- east_model_sum
masked <- Which(trees_proj == 0, cells=TRUE)
tree_data[masked] <- NA

plot(tree_data, main = "Available data for tree areas")
plot(ER_proj, border = 'black', add = TRUE)

tree_data_vec<- as.vector(tree_data)

shrub_data <- east_model_sum
masked <- Which(shrubs_proj == 0, cells=TRUE)
shrub_data[masked] <- NA

plot(shrub_data, main = "Available data for shrub areas")
plot(ER_proj, border = 'black', add = TRUE)

shrub_data_vec<- as.vector(shrub_data)

grass_data <- east_model_sum
masked <- Which(grass_proj == 0, cells=TRUE)
grass_data[masked] <- NA

plot(grass_data, main = "Available data for grass areas")
plot(ER_proj, border = 'black', add = TRUE)

grass_data_vec<- as.vector(grass_data)

clear_data <- east_model_sum
masked <- Which(clear_proj == 0, cells=TRUE)
clear_data[masked] <- NA

plot(clear_data, main = "Available data for clear areas")
plot(ER_proj, border = 'black', add = TRUE)

clear_data_vec<- as.vector(clear_data)

other_data <- east_model_sum
masked <- Which(other_proj == 0, cells=TRUE)
other_data[masked] <- NA

plot(other_data, main = "Available data for other areas")
plot(ER_proj, border = 'black', add = TRUE)

other_data_vec<- as.vector(other_data)


boxplot (other_data_vec, clear_data_vec, grass_data_vec, shrub_data_vec, tree_data_vec,
        main = "Data available by land cover type", 
        at = c(1,2,3,4,5),
        names = c('Other', 'Clear', 'Grass', 'Shrubs', 'Trees'),
        las = 2,
        xlab = "# of Predictions", 
        col = c("#6A3D9A", "#CAB2D6", "#FDBF6F", "#B2DF8A", "#33A02C"),
        horizontal = TRUE)

