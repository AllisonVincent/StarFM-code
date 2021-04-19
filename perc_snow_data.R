setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Analysis_Tests/WY2012/Results')

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


#data_orig<- brick('./2015_03222015_05022015_East_fusion.tif')

#data<- brick('./2015_03222015_05022015_East_fusion.tif')


############ Convert all non-data values for each raster layer to NA, and then convert the values to their binary snow status


# for (i in 1:nlayers(data)) {
#   # Use the raster Which() function for speed:
#   masked <- Which(data[[i]] == -11111, cells=TRUE)
#   data[[i]][masked] <- NA
# }
# 
# 
# for (i in 1:nlayers(data)) {
#   # Use the raster Which() function for speed:
#   masked <- Which(data[[i]] >= 4000, cells=TRUE)
#   data[[i]][masked] <- 1
# }
# 
# 
# for (i in 1:nlayers(data)) {
#   # Use the raster Which() function for speed:
#   masked <- Which(data_orig[[i]] < 4000 & data_orig[[i]] > -11111, cells=TRUE)
#   data[[i]][masked] <- 0
# }


########## Plot the data to confirm the above loops worked correctly

ER<- readOGR("C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/ER_watershed_shp/EastRiver_Project.shp") #shapefile of the study watershed for reference

snow_sum<- raster("./WY2012_snow_sum.tif")

data_proj<- crs(snow_sum)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)


### colors tried: 2012 = brewer YlGnBu , 2010 = topo.colors(150), 2008 = brewer RdBu

dev.new(height=0.91*nrow(snow_sum)/50, width=1.09*ncol(snow_sum)/50)
par(mar = c(5,5,5,3.5))
plot(snow_sum, col = rev(topo.colors(140)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Total Snow Covered Days", cex = 2.0), line = 1.0)
plot(snow_sum, legend.only = TRUE, col = rev(topo.colors(140)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


########## Do some raster math to see where the snow is over time

# snow_sum<- calc(data, fun = sum, na.rm = TRUE)
# 
# plot(snow_sum, main= "Number of snow covered days")
# mtext("22 Mar - 4 May 2015", line = 0.5)
# plot(ER_proj, border = "black", add = TRUE)


######### Compare number of snow-covered days to topographic features

east_DEM<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/East_DEM.tif')

east_slope<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/east_slope.tif')

east_aspect<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_season/DEM_analysis/east_aspect.tif')


######### Elevation statistics

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

levelplot(elev_bands_fac, col.regions = terrain.colors(3), main = "Elevation Bands") + layer(sp.polygons(ER_proj))


#writeRaster(elev_bands, filename = "C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/DEM_Analysis/elev_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)


## Plot the number of correct model instances by elevation band

# For low elevations

# low_elev_snow <- snow_sum
# masked <- Which(elev_low == 0, cells=TRUE)
# low_elev_snow[masked] <- NA
# 
# plot(low_elev_snow, col = brewer.pal(9, name = "YlGnBu"), main = "Snow Days at Low Elevation")
# plot(ER_proj, border = 'black', add = TRUE)

# low_elev_snow_r2<- raster("./low_elev_snow.tif")
# plot(low_elev_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# For mid elevations

# med_elev_snow <- snow_sum
# masked <- Which(elev_med == 0, cells=TRUE)
# med_elev_snow[masked] <- NA
# 
# plot(med_elev_snow, main = "Snow Days at Med Elevation")
# plot(ER_proj, border = 'black', add = TRUE)

# med_elev_snow_r2<- raster("./med_elev_snow.tif")
# plot(med_elev_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# For high elevations

# high_elev_snow <- snow_sum
# masked <- Which(elev_high == 0, cells=TRUE)
# high_elev_snow[masked] <- NA
# 
# plot(high_elev_snow, main = "Snow Days at High Elevation")
# plot(ER_proj, border = 'black', add = TRUE)


# high_elev_snow_r2<- raster("./high_elev_snow.tif")
# plot(high_elev_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)



### Create box and whisker plots to better analyze number of snow days by elevation band

# low_elev_snow_vec<- as.vector(low_elev_snow)
# 
# mid_elev_snow_vec<- as.vector(med_elev_snow)
# 
# high_elev_snow_vec<- as.vector(high_elev_snow)
# 
# 
# boxplot(low_elev_snow_vec, mid_elev_snow_vec, high_elev_snow_vec, 
#         main = "Snow-covered days at elevation bands", 
#         at = c(1,2,3),
#         names = c('Low', 'Mid', 'High'),
#         las = 2,
#         xlab = "# of snow days", 
#         col = brewer.pal(n = 3, name = "Dark2"), 
#         horizontal = TRUE)



############# Analyses by slope

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

levelplot(slope_bands_fac, col.regions = rev(brewer.pal(3, name = "RdYlBu")), main = "Slope Bands") + layer(sp.polygons(ER_proj))


#writeRaster(slope_bands, filename = "C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/DEM_Analysis/slope_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)

## Plot the number of correct model instances by slope band

# For low slopes

# low_slope_snow <- snow_sum
# masked <- Which(slope_low == 0, cells=TRUE)
# low_slope_snow[masked] <- NA
# 
# low_slope_snow_crop<- crop(low_slope_snow, extent(low_slope_snow, 2, 1590, 2, 1805))
# 
# plot(low_slope_snow_crop, main = "Snow Days at Low Slopes")
# plot(ER_proj, border = 'black', add = TRUE)

# low_slope_snow_r2<- raster("./low_slope_snow.tif")
# plot(low_slope_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)



# For mid slopes

# med_slope_snow <- snow_sum
# masked <- Which(slope_med == 0, cells=TRUE)
# med_slope_snow[masked] <- NA
# 
# med_slope_snow_crop<- crop(med_slope_snow, extent(med_slope_snow, 2, 1590, 2, 1805))
# 
# plot(med_slope_snow_crop, main = "Snow Days at Med Slopes")
# plot(ER_proj, border = 'black', add = TRUE)

# med_slope_snow_r2<- raster("./med_slope_snow.tif")
# plot(med_slope_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)



# For high slopes

# high_slope_snow <- snow_sum
# masked <- Which(slope_high == 0, cells=TRUE)
# high_slope_snow[masked] <- NA
# 
# high_slope_snow_crop<- crop(high_slope_snow, extent(high_slope_snow, 2, 1590, 2, 1805))
# 
# plot(high_slope_snow_crop, main = "Snow Days at High Slopes")
# plot(ER_proj, border = 'black', add = TRUE)


# high_slope_snow_r2<- raster("./high_slope_snow.tif")
# plot(high_slope_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


### Create box and whisker plots to better analyze number of snow days by slope band

# low_slope_snow_vec<- as.vector(low_slope_snow_crop)
# 
# mid_slope_snow_vec<- as.vector(med_slope_snow_crop)
# 
# high_slope_snow_vec<- as.vector(high_slope_snow_crop)
# 
# 
# boxplot(low_slope_snow_vec, mid_slope_snow_vec, high_slope_snow_vec, 
#         main = "Snow-covered days at slope bands", 
#         at = c(1,2,3),
#         names = c('Low', 'Mid', 'High'),
#         las = 2,
#         xlab = "# of snow days", 
#         col = brewer.pal(n = 3, name = "Dark2"), 
#         horizontal = TRUE)


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

levelplot(aspect_cat, att = 'direction', col.regions = c("#FF0000", "#FFA500", "#F0E68C", "#87CEEB"), main = 'Aspect') + layer(sp.polygons(ER_proj))

#writeRaster(aspect_cat, filename = "C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/DEM_Analysis/aspect_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)

## Plot the number of correct model instances by direction

# For north facing slopes

# north_aspect_snow <- snow_sum
# masked <- Which(north == 0, cells=TRUE)
# north_aspect_snow[masked] <- NA
# 
# north_aspect_snow_crop<- crop(north_aspect_snow, extent(north_aspect_snow, 2, 1590, 2, 1805))
# 
# plot(north_aspect_snow_crop, main = "Snow days at north facing slopes")
# plot(ER_proj, border = 'black', add = TRUE)


# north_aspect_snow_r2<- raster("./north_aspect_snow.tif")
# plot(north_aspect_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# For east facing slopes

# east_aspect_snow <- snow_sum
# masked <- Which(east == 0, cells=TRUE)
# east_aspect_snow[masked] <- NA
# 
# east_aspect_snow_crop<- crop(east_aspect_snow, extent(east_aspect_snow, 2, 1590, 2, 1805))
# 
# plot(east_aspect_snow_crop, main = "Snow days at east facing slopes")
# plot(ER_proj, border = 'black', add = TRUE)


# east_aspect_snow_r2<- raster("./east_aspect_snow.tif")
# plot(east_aspect_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# For south facing slopes

# south_aspect_snow <- snow_sum
# masked <- Which(south == 0, cells=TRUE)
# south_aspect_snow[masked] <- NA
# 
# south_aspect_snow_crop<- crop(south_aspect_snow, extent(south_aspect_snow, 2, 1590, 2, 1805))
# 
# plot(south_aspect_snow_crop, main = "Snow days at south facing slopes")
# plot(ER_proj, border = 'black', add = TRUE)


# south_aspect_snow_r2<- raster("./south_aspect_snow.tif")
# plot(south_aspect_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# For west facing slopes

# west_aspect_snow <- snow_sum
# masked <- Which(west == 0, cells=TRUE)
# west_aspect_snow[masked] <- NA
# 
# west_aspect_snow_crop<- crop(west_aspect_snow, extent(west_aspect_snow, 2, 1590, 2, 1805))
# 
# plot(west_aspect_snow_crop, main = "Snow days at west facing slopes")
# plot(ER_proj, border = 'black', add = TRUE)


# west_aspect_snow_r2<- raster("./west_aspect_snow.tif")
# plot(west_aspect_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


### Create box and whisker plots to better analyze data availability by aspect direction


# north_aspect_snow_vec<- as.vector(north_aspect_snow_crop)
# 
# east_aspect_snow_vec<- as.vector(east_aspect_snow_crop)
# 
# south_aspect_snow_vec<- as.vector(south_aspect_snow_crop)
# 
# west_aspect_snow_vec<- as.vector(west_aspect_snow_crop)
# 
# 
# boxplot(north_aspect_snow_vec, east_aspect_snow_vec, south_aspect_snow_vec, west_aspect_snow_vec, 
#         main = "Snow-covered days at aspect", 
#         at = c(1,2,3,4),
#         names = c('North', 'East', 'South', 'West'),
#         las = 2,
#         xlab = "# of snow days", 
#         col = brewer.pal(n = 4, name = "Dark2"), 
#         horizontal = TRUE)


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


east_trees <- (east_veg_class == 1)
east_shrub <- (east_veg_class == 2)
east_grass <- (east_veg_class == 3 | east_veg_class == 8)
east_clear <- (east_veg_class == 5 | east_veg_class == 9)
east_other<- (east_veg_class == 4 | east_veg_class == 6 | east_veg_class == 7)


## Same as above, just all in one raster
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

#writeRaster(veg_class, filename = "C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Validation_Tests/Full_Season/DEM_Analysis/landcover_bands.tif", bandorder='BSQ', datatype='INT2S', formatt='GTiff', overwrite=TRUE)

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


# trees_snow <- snow_sum
# masked <- Which(trees_proj == 0, cells = TRUE)
# trees_snow[masked] <- NA
# 
# plot(trees_snow, main = "Snow for Tree Areas")
# plot(ER_proj, border = 'black', add = TRUE)


# trees_snow_r2<- raster("./trees_snow.tif")
# plot(trees_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)



# shrubs_snow <- snow_sum
# masked <- Which(shrubs_proj == 0, cells = TRUE)
# shrubs_snow[masked] <- NA
# 
# plot(shrubs_snow, main = "Snow for Shrub Areas")
# plot(ER_proj, border = 'black', add = TRUE)

# shrubs_snow_r2<- raster("./shrubs_snow.tif")
# plot(shrubs_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)



# grass_snow <- snow_sum
# masked <- Which(grass_proj == 0, cells = TRUE)
# grass_snow[masked] <- NA
# 
# plot(grass_snow, main = "Snow for Grass Areas")
# plot(ER_proj, border = 'black', add = TRUE)

# grass_snow_r2<- raster("./grass_snow.tif")
# plot(grass_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# clear_snow <- snow_sum
# masked <- Which(clear_proj == 0, cells=TRUE)
# clear_snow[masked] <- NA
# 
# plot(clear_snow, main = "Snow for Clear Areas")
# plot(ER_proj, border = 'black', add = TRUE)

# clear_snow_r2<- raster("./clear_snow.tif")
# plot(clear_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)


# other_snow <- snow_sum
# masked <- Which(other_proj == 0, cells=TRUE)
# other_snow[masked] <- NA
# 
# plot(other_snow, main = "Snow for All Other Areas")
# plot(ER_proj, border = 'black', add = TRUE)

# other_snow_r2<- raster("./other_snow.tif")
# plot(other_snow_r2)
# plot(ER_proj, border = 'black', add = TRUE)

######### Create a boxplot for snow 

# trees_snow_vec<- as.vector(trees_snow)
# 
# shrubs_snow_vec<- as.vector(shrubs_snow)
# 
# grass_snow_vec<- as.vector(grass_snow)
# 
# clear_snow_vec<- as.vector(clear_snow)
# 
# other_snow_vec<- as.vector(other_snow)
# 
# 
# boxplot(trees_snow_vec, shrubs_snow_vec, grass_snow_vec, clear_snow_vec, other_snow_vec, 
#         main = "Snow-covered days for landcover class", 
#         at = c(1,2,3,4,5),
#         names = c('trees', 'shrubs', 'grass', 'clear', 'other'),
#         las = 2,
#         xlab = "# of snow days", 
#         col = brewer.pal(n = 5, name = "Dark2"), 
#         horizontal = TRUE)


#################### Percentage of snow covered days

setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Analysis_Tests/WY2008/Results')

total_model<- raster('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Analysis_Tests/WY2008/WY2008_model_sum.tif')

snow_sum<- raster("./WY2008_snow_sum.tif")



plot(total_model, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Total Snow Covered Days", cex = 2.0), line = 1.0)
plot(snow_sum, legend.only = TRUE, col = rev(topo.colors(140)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


corrected<- (snow_sum/total_model) * 100


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(corrected, col = rev(topo.colors(100)),  main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Annual Snow Covered Days", cex = 1.5), line = 1.0)
plot(corrected, legend.only = TRUE, col = rev(topo.colors(140)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



##### By elevation

# For low elevations

perc_low_elev_snow <- corrected
masked <- Which(elev_low == 0, cells=TRUE)
perc_low_elev_snow[masked] <- NA


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_low_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Low Elev", cex = 1.5), line = 1.0)
plot(perc_low_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For mid elevations

perc_med_elev_snow <- corrected
masked <- Which(elev_med == 0, cells=TRUE)
perc_med_elev_snow[masked] <- NA


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_med_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Mid Elev", cex = 1.5), line = 1.0)
plot(perc_med_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For high elevations

perc_high_elev_snow <- corrected
masked <- Which(elev_high == 0, cells=TRUE)
perc_high_elev_snow[masked] <- NA


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_high_elev_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at High Elev", cex = 1.5), line = 1.0)
plot(perc_high_elev_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


## boxplot for percent snow covered days by elevation

low_elev_snow_vec<- as.vector(perc_low_elev_snow)

mid_elev_snow_vec<- as.vector(perc_med_elev_snow)

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
masked <- Which(slope_low == 0, cells=TRUE)
perc_low_slope_snow[masked] <- NA

perc_low_slope_snow_crop<-crop(perc_low_slope_snow, extent(perc_low_slope_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_low_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Low Slope", cex = 1.5), line = 1.0)
plot(perc_low_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

# For mid slopes

perc_med_slope_snow <- corrected
masked <- Which(slope_med == 0, cells=TRUE)
perc_med_slope_snow[masked] <- NA

perc_med_slope_snow_crop<- crop(perc_med_slope_snow, extent(perc_med_slope_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_med_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at Mid Slope", cex = 1.5), line = 1.0)
plot(perc_med_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

# For high slopes

perc_high_slope_snow <- corrected
masked <- Which(slope_high == 0, cells=TRUE)
perc_high_slope_snow[masked] <- NA

perc_high_slope_snow_crop<- crop(perc_high_slope_snow, extent(perc_high_slope_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_high_slope_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at High Slope", cex = 1.5), line = 1.0)
plot(perc_high_slope_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


### Create box and whisker plots to better analyze number of snow days by elevation band

low_slope_snow_vec<- as.vector(perc_low_slope_snow_crop)

mid_slope_snow_vec<- as.vector(perc_med_slope_snow_crop)

high_slope_snow_vec<- as.vector(perc_high_slope_snow_crop)


boxplot(low_slope_snow_vec, mid_slope_snow_vec, high_slope_snow_vec, 
        main = "WY 2008 Percent Snow-Covered Days at Slope Band", cex.main = 1.2,
        at = c(1,2,3),
        names = c('Low', 'Mid', 'High'), cex.names = 1.2, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = "#7570B3", 
        horizontal = TRUE)

############ By aspect


# For north facing slopes

perc_north_aspect_snow <- corrected
masked <- Which(north == 0, cells=TRUE)
perc_north_aspect_snow[masked] <- NA

perc_north_aspect_snow_crop<- crop(perc_north_aspect_snow, extent(perc_north_aspect_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_north_aspect_snow_crop, col = rev(topo.colors(150)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at North Aspect", cex = 1.5), line = 1.0)
plot(perc_north_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For east facing slopes

perc_east_aspect_snow <- corrected
masked <- Which(east == 0, cells=TRUE)
perc_east_aspect_snow[masked] <- NA

perc_east_aspect_snow_crop<- crop(perc_east_aspect_snow, extent(perc_east_aspect_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_east_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at East Aspect", cex = 1.5), line = 1.0)
plot(perc_east_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For south facing slopes

perc_south_aspect_snow <- corrected
masked <- Which(south == 0, cells=TRUE)
perc_south_aspect_snow[masked] <- NA

perc_south_aspect_snow_crop<- crop(perc_south_aspect_snow, extent(perc_south_aspect_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_south_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at South Aspect", cex = 1.5), line = 1.0)
plot(perc_south_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


# For west facing slopes

perc_west_aspect_snow <- corrected
masked <- Which(west == 0, cells=TRUE)
perc_west_aspect_snow[masked] <- NA

perc_west_aspect_snow_crop<- crop(perc_west_aspect_snow, extent(perc_west_aspect_snow, 2, 1590, 2, 1805))


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_west_aspect_snow_crop, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days at West Aspect", cex = 1.5), line = 1.0)
plot(perc_west_aspect_snow_crop, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)

### Create box and whisker plots to better analyze data availability by aspect direction


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

perc_trees_snow <- corrected
masked <- Which(trees_proj == 0, cells = TRUE)
perc_trees_snow[masked] <- NA


dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_trees_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Tree Areas", cex = 1.5), line = 1.0)
plot(perc_trees_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



perc_shrubs_snow <- corrected
masked <- Which(shrubs_proj == 0, cells = TRUE)
perc_shrubs_snow[masked] <- NA

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_shrubs_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Shrub Areas", cex = 1.5), line = 1.0)
plot(perc_shrubs_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)



perc_grass_snow <- corrected
masked <- Which(grass_proj == 0, cells = TRUE)
perc_grass_snow[masked] <- NA

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_grass_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Grass Areas", cex = 1.5), line = 1.0)
plot(perc_grass_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


perc_clear_snow <- corrected
masked <- Which(clear_proj == 0, cells=TRUE)
perc_clear_snow[masked] <- NA

dev.new(height=0.91*nrow(corrected)/50, width=1.09*ncol(corrected)/50)
par(mar = c(5,5,5,3.5))
plot(perc_clear_snow, col = rev(topo.colors(100)), main = "", xlab = "Longitude", ylab = "Latitude", cex.lab = 2, cex.axis = 2, legend = FALSE)
title(main = list("WY 2012 Percent Snow Days for Clear Areas", cex = 1.5), line = 1.0)
plot(perc_clear_snow, legend.only = TRUE, col = rev(topo.colors(100)), axis.args = list(cex.axis = 1.9))
plot(ER_proj, border = "black", lwd = 2, add = TRUE)


perc_other_snow <- corrected
masked <- Which(other_proj == 0, cells=TRUE)
perc_other_snow[masked] <- NA

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_other_snow, col = rev(topo.colors(100)), main = "WY 2010 percent snow for all other areas")
plot(ER_proj, border = 'black', add = TRUE)


######### Create a boxplot for snow 

trees_snow_vec<- as.vector(perc_trees_snow)

shrubs_snow_vec<- as.vector(perc_shrubs_snow)

grass_snow_vec<- as.vector(perc_grass_snow)

clear_snow_vec<- as.vector(perc_clear_snow)

other_snow_vec<- as.vector(perc_other_snow)


boxplot(trees_snow_vec, shrubs_snow_vec, grass_snow_vec, clear_snow_vec, 
        main = "WY 2012 Percent Snow-Covered Days by Landcover Class", cex.main = 1.2,
        at = c(1,2,3,4),
        names = c('Trees', 'Shrubs', 'Grass', 'Clear'), cex.names = 1.2, cex.axis = 1.1,
        las = 2,
        xlab = "percent snow days", cex.lab = 1.3,
        col = "#D95F02", 
        horizontal = TRUE)



###### Create boxplots for slope by elevation band to isolate the elevation relationship between snow cover and landscape characteristics

## low slopes/low elevation


perc_low_slope_low_elev<- (elev_low == 1 & slope_low == 1)
perc_low_slope_low_elev_vec<- as.vector(perc_low_slope_low_elev)
perc_low_slope_low_elev_vec_sum<- sum(perc_low_slope_low_elev_vec, na.rm = TRUE)

perc_low_slope_low_elev_snow <- corrected
masked <- Which(perc_low_slope_low_elev == 0, cells=TRUE)
perc_low_slope_low_elev_snow[masked] <- NA

perc_low_slope_low_elev_snow_crop<- crop(perc_low_slope_low_elev_snow, extent(perc_low_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2008 percent snow days at low slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## low slopes/mid elevation


perc_low_slope_mid_elev<- (elev_med == 1 & slope_low == 1)
perc_low_slope_mid_elev_vec<- as.vector(perc_low_slope_mid_elev)
perc_low_slope_mid_elev_vec_sum<- sum(perc_low_slope_mid_elev_vec, na.rm = TRUE)

perc_low_slope_mid_elev_snow <- corrected
masked <- Which(perc_low_slope_mid_elev == 0, cells=TRUE)
perc_low_slope_mid_elev_snow[masked] <- NA

perc_low_slope_mid_elev_snow_crop<- crop(perc_low_slope_mid_elev_snow, extent(perc_low_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2008 percent snow days at low slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## low slopes/high elevation


perc_low_slope_high_elev<- (elev_high == 1 & slope_low == 1)
perc_low_slope_high_elev_vec<- as.vector(perc_low_slope_high_elev)
perc_low_slope_high_elev_vec_sum<- sum(perc_low_slope_high_elev_vec, na.rm = TRUE)

perc_low_slope_high_elev_snow <- corrected
masked <- Which(perc_low_slope_high_elev == 0, cells=TRUE)
perc_low_slope_high_elev_snow[masked] <- NA

perc_low_slope_high_elev_snow_crop<- crop(perc_low_slope_high_elev_snow, extent(perc_low_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_low_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2008 percent snow days at low slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    



### create boxplots of the above

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


perc_mid_slope_low_elev<- (elev_low == 1 & slope_med == 1)
perc_mid_slope_low_elev_vec<- as.vector(perc_mid_slope_low_elev)
perc_mid_slope_low_elev_vec_sum<- sum(perc_mid_slope_low_elev_vec, na.rm = TRUE)

perc_mid_slope_low_elev_snow <- corrected
masked <- Which(perc_mid_slope_low_elev == 0, cells=TRUE)
perc_mid_slope_low_elev_snow[masked] <- NA

perc_mid_slope_low_elev_snow_crop<- crop(perc_mid_slope_low_elev_snow, extent(perc_mid_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## mid slopes/mid elevation


perc_mid_slope_mid_elev<- (elev_med == 1 & slope_med == 1)
perc_mid_slope_mid_elev_vec<- as.vector(perc_mid_slope_mid_elev)
perc_mid_slope_mid_elev_vec_sum<- sum(perc_mid_slope_mid_elev_vec, na.rm = TRUE)

perc_mid_slope_mid_elev_snow <- corrected
masked <- Which(perc_mid_slope_mid_elev == 0, cells=TRUE)
perc_mid_slope_mid_elev_snow[masked] <- NA

perc_mid_slope_mid_elev_snow_crop<- crop(perc_mid_slope_mid_elev_snow, extent(perc_mid_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## mid slopes/high elevation


perc_mid_slope_high_elev<- (elev_high == 1 & slope_med == 1)
perc_mid_slope_high_elev_vec<- as.vector(perc_mid_slope_high_elev)
perc_mid_slope_high_elev_vec_sum<- sum(perc_mid_slope_high_elev_vec, na.rm = TRUE)

perc_mid_slope_high_elev_snow <- corrected
masked <- Which(perc_mid_slope_high_elev == 0, cells=TRUE)
perc_mid_slope_high_elev_snow[masked] <- NA

perc_mid_slope_high_elev_snow_crop<- crop(perc_mid_slope_high_elev_snow, extent(perc_mid_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_mid_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at mid slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    



### create boxplots of the above

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


perc_high_slope_low_elev<- (elev_low == 1 & slope_high == 1)
perc_high_slope_low_elev_vec<- as.vector(perc_high_slope_low_elev)
perc_high_slope_low_elev_vec_sum<- sum(perc_high_slope_low_elev_vec, na.rm = TRUE)

perc_high_slope_low_elev_snow <- corrected
masked <- Which(perc_high_slope_low_elev == 0, cells=TRUE)
perc_high_slope_low_elev_snow[masked] <- NA

perc_high_slope_low_elev_snow_crop<- crop(perc_high_slope_low_elev_snow, extent(perc_high_slope_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## high slopes/mid elevation


perc_high_slope_mid_elev<- (elev_med == 1 & slope_high == 1)
perc_high_slope_mid_elev_vec<- as.vector(perc_high_slope_mid_elev)
perc_high_slope_mid_elev_vec_sum<- sum(perc_high_slope_mid_elev_vec, na.rm = TRUE)

perc_high_slope_mid_elev_snow <- corrected
masked <- Which(perc_high_slope_mid_elev == 0, cells=TRUE)
perc_high_slope_mid_elev_snow[masked] <- NA

perc_high_slope_mid_elev_snow_crop<- crop(perc_high_slope_mid_elev_snow, extent(perc_high_slope_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)  


## high slopes/high elevation


perc_high_slope_high_elev<- (elev_high == 1 & slope_high == 1)
perc_high_slope_high_elev_vec<- as.vector(perc_high_slope_high_elev)
perc_high_slope_high_elev_vec_sum<- sum(perc_high_slope_high_elev_vec, na.rm = TRUE)

perc_high_slope_high_elev_snow <- corrected
masked <- Which(perc_high_slope_high_elev == 0, cells=TRUE)
perc_high_slope_high_elev_snow[masked] <- NA

perc_high_slope_high_elev_snow_crop<- crop(perc_high_slope_high_elev_snow, extent(perc_high_slope_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_high_slope_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at high slope/high elevations")
plot(ER_proj, border = 'black', add = TRUE)    



### create boxplots of the above

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

perc_north_aspect_low_elev<- (elev_low == 1 & north == 1)
perc_north_aspect_low_elev_vec<- as.vector(perc_north_aspect_low_elev)
perc_north_aspect_low_elev_vec_sum<- sum(perc_north_aspect_low_elev_vec, na.rm = TRUE)

perc_north_aspect_low_elev_snow <- corrected
masked <- Which(perc_north_aspect_low_elev == 0, cells=TRUE)
perc_north_aspect_low_elev_snow[masked] <- NA

perc_north_aspect_low_elev_snow_crop<- crop(perc_north_aspect_low_elev_snow, extent(perc_north_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## north aspect/mid elevation

perc_north_aspect_mid_elev<- (elev_med == 1 & north == 1)
perc_north_aspect_mid_elev_vec<- as.vector(perc_north_aspect_mid_elev)
perc_north_aspect_mid_elev_vec_sum<- sum(perc_north_aspect_mid_elev_vec, na.rm = TRUE)

perc_north_aspect_mid_elev_snow <- corrected
masked <- Which(perc_north_aspect_mid_elev == 0, cells=TRUE)
perc_north_aspect_mid_elev_snow[masked] <- NA

perc_north_aspect_mid_elev_snow_crop<- crop(perc_north_aspect_mid_elev_snow, extent(perc_north_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## north aspect/high elevation

perc_north_aspect_high_elev<- (elev_high == 1 & north == 1)
perc_north_aspect_high_elev_vec<- as.vector(perc_north_aspect_high_elev)
perc_north_aspect_high_elev_vec_sum<- sum(perc_north_aspect_high_elev_vec, na.rm = TRUE)

perc_north_aspect_high_elev_snow <- corrected
masked <- Which(perc_north_aspect_high_elev == 0, cells=TRUE)
perc_north_aspect_high_elev_snow[masked] <- NA

perc_north_aspect_high_elev_snow_crop<- crop(perc_north_aspect_high_elev_snow, extent(perc_north_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_north_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at north aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_east_aspect_low_elev<- (elev_low == 1 & east == 1)
perc_east_aspect_low_elev_vec<- as.vector(perc_east_aspect_low_elev)
perc_east_aspect_low_elev_vec_sum<- sum(perc_east_aspect_low_elev_vec, na.rm = TRUE)

perc_east_aspect_low_elev_snow <- corrected
masked <- Which(perc_east_aspect_low_elev == 0, cells=TRUE)
perc_east_aspect_low_elev_snow[masked] <- NA

perc_east_aspect_low_elev_snow_crop<- crop(perc_east_aspect_low_elev_snow, extent(perc_east_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## east aspect/mid elevation

perc_east_aspect_mid_elev<- (elev_med == 1 & east == 1)
perc_east_aspect_mid_elev_vec<- as.vector(perc_east_aspect_mid_elev)
perc_east_aspect_mid_elev_vec_sum<- sum(perc_east_aspect_mid_elev_vec, na.rm = TRUE)

perc_east_aspect_mid_elev_snow <- corrected
masked <- Which(perc_east_aspect_mid_elev == 0, cells=TRUE)
perc_east_aspect_mid_elev_snow[masked] <- NA

perc_east_aspect_mid_elev_snow_crop<- crop(perc_east_aspect_mid_elev_snow, extent(perc_east_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## east aspect/high elevation

perc_east_aspect_high_elev<- (elev_high == 1 & east == 1)
perc_east_aspect_high_elev_vec<- as.vector(perc_east_aspect_high_elev)
perc_east_aspect_high_elev_vec_sum<- sum(perc_east_aspect_high_elev_vec, na.rm = TRUE)

perc_east_aspect_high_elev_snow <- corrected
masked <- Which(perc_east_aspect_high_elev == 0, cells=TRUE)
perc_east_aspect_high_elev_snow[masked] <- NA

perc_east_aspect_high_elev_snow_crop<- crop(perc_east_aspect_high_elev_snow, extent(perc_east_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_east_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at east aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_south_aspect_low_elev<- (elev_low == 1 & south == 1)
perc_south_aspect_low_elev_vec<- as.vector(perc_south_aspect_low_elev)
perc_south_aspect_low_elev_vec_sum<- sum(perc_south_aspect_low_elev_vec, na.rm = TRUE)

perc_south_aspect_low_elev_snow <- corrected
masked <- Which(perc_south_aspect_low_elev == 0, cells=TRUE)
perc_south_aspect_low_elev_snow[masked] <- NA

perc_south_aspect_low_elev_snow_crop<- crop(perc_south_aspect_low_elev_snow, extent(perc_south_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## south aspect/mid elevation

perc_south_aspect_mid_elev<- (elev_med == 1 & south == 1)
perc_south_aspect_mid_elev_vec<- as.vector(perc_south_aspect_mid_elev)
perc_south_aspect_mid_elev_vec_sum<- sum(perc_south_aspect_mid_elev_vec, na.rm = TRUE)

perc_south_aspect_mid_elev_snow <- corrected
masked <- Which(perc_south_aspect_mid_elev == 0, cells=TRUE)
perc_south_aspect_mid_elev_snow[masked] <- NA

perc_south_aspect_mid_elev_snow_crop<- crop(perc_south_aspect_mid_elev_snow, extent(perc_south_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## south aspect/high elevation

perc_south_aspect_high_elev<- (elev_high == 1 & south == 1)
perc_south_aspect_high_elev_vec<- as.vector(perc_south_aspect_high_elev)
perc_south_aspect_high_elev_vec_sum<- sum(perc_south_aspect_high_elev_vec, na.rm = TRUE)

perc_south_aspect_high_elev_snow <- corrected
masked <- Which(perc_south_aspect_high_elev == 0, cells=TRUE)
perc_south_aspect_high_elev_snow[masked] <- NA

perc_south_aspect_high_elev_snow_crop<- crop(perc_south_aspect_high_elev_snow, extent(perc_south_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_south_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at south aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_west_aspect_low_elev<- (elev_low == 1 & west == 1)
perc_west_aspect_low_elev_vec<- as.vector(perc_west_aspect_low_elev)
perc_west_aspect_low_elev_vec_sum<- sum(perc_west_aspect_low_elev_vec, na.rm = TRUE)

perc_west_aspect_low_elev_snow <- corrected
masked <- Which(perc_west_aspect_low_elev == 0, cells=TRUE)
perc_west_aspect_low_elev_snow[masked] <- NA

perc_west_aspect_low_elev_snow_crop<- crop(perc_west_aspect_low_elev_snow, extent(perc_west_aspect_low_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_low_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## west aspect/mid elevation

perc_west_aspect_mid_elev<- (elev_med == 1 & west == 1)
perc_west_aspect_mid_elev_vec<- as.vector(perc_west_aspect_mid_elev)
perc_west_aspect_mid_elev_vec_sum<- sum(perc_west_aspect_mid_elev_vec, na.rm = TRUE)

perc_west_aspect_mid_elev_snow <- corrected
masked <- Which(perc_west_aspect_mid_elev == 0, cells=TRUE)
perc_west_aspect_mid_elev_snow[masked] <- NA

perc_west_aspect_mid_elev_snow_crop<- crop(perc_west_aspect_mid_elev_snow, extent(perc_west_aspect_mid_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_mid_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## west aspect/high elevation

perc_west_aspect_high_elev<- (elev_high == 1 & west == 1)
perc_west_aspect_high_elev_vec<- as.vector(perc_west_aspect_high_elev)
perc_west_aspect_high_elev_vec_sum<- sum(perc_west_aspect_high_elev_vec, na.rm = TRUE)

perc_west_aspect_high_elev_snow <- corrected
masked <- Which(perc_west_aspect_high_elev == 0, cells=TRUE)
perc_west_aspect_high_elev_snow[masked] <- NA

perc_west_aspect_high_elev_snow_crop<- crop(perc_west_aspect_high_elev_snow, extent(perc_west_aspect_high_elev_snow, 2, 1590, 2, 1805))

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_west_aspect_high_elev_snow_crop, col = rev(topo.colors(100)), main = "WY 2012 percent snow days at west aspect/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_trees_low_elev<- (elev_low == 1 & trees_proj == 1)
perc_trees_low_elev_vec<- as.vector(perc_trees_low_elev)
perc_trees_low_elev_vec_sum<- sum(perc_trees_low_elev_vec, na.rm = TRUE)

perc_trees_low_elev_snow <- corrected
masked <- Which(perc_trees_low_elev == 0, cells=TRUE)
perc_trees_low_elev_snow[masked] <- NA

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## tree areas/mid elevation

perc_trees_mid_elev<- (elev_med == 1 & trees_proj == 1)
perc_trees_mid_elev_vec<- as.vector(perc_trees_mid_elev)
perc_trees_mid_elev_vec_sum<- sum(perc_trees_mid_elev_vec, na.rm = TRUE)

perc_trees_mid_elev_snow <- corrected
masked <- Which(perc_trees_mid_elev == 0, cells=TRUE)
perc_trees_mid_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## tree areas/high elevation

perc_trees_high_elev<- (elev_high == 1 & trees_proj == 1)
perc_trees_high_elev_vec<- as.vector(perc_trees_high_elev)
perc_trees_high_elev_vec_sum<- sum(perc_trees_high_elev_vec, na.rm = TRUE)

perc_trees_high_elev_snow <- corrected
masked <- Which(perc_trees_high_elev == 0, cells=TRUE)
perc_trees_high_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_trees_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for tree areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_shrubs_low_elev<- (elev_low == 1 & shrubs_proj == 1)
perc_shrubs_low_elev_vec<- as.vector(perc_shrubs_low_elev)
perc_shrubs_low_elev_vec_sum<- sum(perc_shrubs_low_elev_vec, na.rm = TRUE)

perc_shrubs_low_elev_snow <- corrected
masked <- Which(perc_shrubs_low_elev == 0, cells=TRUE)
perc_shrubs_low_elev_snow[masked] <- NA

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## shrub areas/mid elevation

perc_shrubs_mid_elev<- (elev_med == 1 & shrubs_proj == 1)
perc_shrubs_mid_elev_vec<- as.vector(perc_shrubs_mid_elev)
perc_shrubs_mid_elev_vec_sum<- sum(perc_shrubs_mid_elev_vec, na.rm = TRUE)

perc_shrubs_mid_elev_snow <- corrected
masked <- Which(perc_shrubs_mid_elev == 0, cells=TRUE)
perc_shrubs_mid_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## shrub areas/high elevation

perc_shrubs_high_elev<- (elev_high == 1 & shrubs_proj == 1)
perc_shrubs_high_elev_vec<- as.vector(perc_shrubs_high_elev)
perc_shrubs_high_elev_vec_sum<- sum(perc_shrubs_high_elev_vec, na.rm = TRUE)

perc_shrubs_high_elev_snow <- corrected
masked <- Which(perc_shrubs_high_elev == 0, cells=TRUE)
perc_shrubs_high_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_shrubs_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for shrub areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_grass_low_elev<- (elev_low == 1 & grass_proj == 1)
perc_grass_low_elev_vec<- as.vector(perc_grass_low_elev)
perc_grass_low_elev_vec_sum<- sum(perc_grass_low_elev_vec, na.rm = TRUE)

perc_grass_low_elev_snow <- corrected
masked <- Which(perc_grass_low_elev == 0, cells=TRUE)
perc_grass_low_elev_snow[masked] <- NA

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## grass areas/mid elevation

perc_grass_mid_elev<- (elev_med == 1 & grass_proj == 1)
perc_grass_mid_elev_vec<- as.vector(perc_grass_mid_elev)
perc_grass_mid_elev_vec_sum<- sum(perc_grass_mid_elev_vec, na.rm = TRUE)

perc_grass_mid_elev_snow <- corrected
masked <- Which(perc_grass_mid_elev == 0, cells=TRUE)
perc_grass_mid_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## grass areas/high elevation

perc_grass_high_elev<- (elev_high == 1 & grass_proj == 1)
perc_grass_high_elev_vec<- as.vector(perc_grass_high_elev)
perc_grass_high_elev_vec_sum<- sum(perc_grass_high_elev_vec, na.rm = TRUE)

perc_grass_high_elev_snow <- corrected
masked <- Which(perc_grass_high_elev == 0, cells=TRUE)
perc_grass_high_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_grass_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for grass areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

perc_clear_low_elev<- (elev_low == 1 & clear_proj == 1)
perc_clear_low_elev_vec<- as.vector(perc_clear_low_elev)
perc_clear_low_elev_vec_sum<- sum(perc_clear_low_elev_vec, na.rm = TRUE)

perc_clear_low_elev_snow <- corrected
masked <- Which(perc_clear_low_elev == 0, cells=TRUE)
perc_clear_low_elev_snow[masked] <- NA

par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_low_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/low elevations")
plot(ER_proj, border = 'black', add = TRUE)


## clear areas/mid elevation

perc_clear_mid_elev<- (elev_med == 1 & clear_proj == 1)
perc_clear_mid_elev_vec<- as.vector(perc_clear_mid_elev)
perc_clear_mid_elev_vec_sum<- sum(perc_clear_mid_elev_vec, na.rm = TRUE)

perc_clear_mid_elev_snow <- corrected
masked <- Which(perc_clear_mid_elev == 0, cells=TRUE)
perc_clear_mid_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_mid_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/mid elevations")
plot(ER_proj, border = 'black', add = TRUE)


## tree areas/high elevation

perc_clear_high_elev<- (elev_high == 1 & clear_proj == 1)
perc_clear_high_elev_vec<- as.vector(perc_clear_high_elev)
perc_clear_high_elev_vec_sum<- sum(perc_clear_high_elev_vec, na.rm = TRUE)

perc_clear_high_elev_snow <- corrected
masked <- Which(perc_clear_high_elev == 0, cells=TRUE)
perc_clear_high_elev_snow[masked] <- NA


par(mar = c(2.5,2.5,2.5,2.5))
plot(perc_clear_high_elev_snow, col = rev(topo.colors(100)), main = "WY 2012 percent snow days for clear areas/high elevations")
plot(ER_proj, border = 'black', add = TRUE)


### create boxplots of the above

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

