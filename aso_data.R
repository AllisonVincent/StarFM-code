### This script is for plotting ASO snow depth data and using it for validation analysis on STARFM data via a confusion matrix

### ASO data can be downloaded from the NASA Earth Data archive (earthdata.nasa.gov)


### Libraries needed
library(sp)
library(sf)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(caret)
library(data.table)
library(rasterVis)



### Load in the ASO data and watershed shapefile

ER<- readOGR('./EastRiver_Project.shp') #shapefile of the study watershed for reference
aso_50m<- raster('./ASO_50M.tif')
aso_3m<- raster('./ASO_3M.tif')

### Check the projections of the data

proj4string(aso_50m)
proj4string(aso_3m)


### Project the ER shapefile to same projection as aso data
data_proj<- crs(aso_50m)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

## Plot the 50-m ASO data as is
#dev.new()
par(mar = c(4,4,4,3.5))
plot(aso_50m, col = rev(viridis(25)), main = "", xlab = "Meters", ylab = "Meters", cex.lab = 1, cex.axis = 0.9, legend = FALSE)
title(main = list("ASO 50-m Snow Depths", cex = 1.5), line = 1.5)
mtext("Date", line = 0.15, cex = 1.2)
plot(aso_50m, legend.only = TRUE, col = rev(viridis(25)), axis.args = list(cex.axis = 1.1), legend.args = list(text = 'Depths (m)', line = 0.5))
plot(ER_proj, border = "black", add = TRUE)



## Plot the 3-m ASO data as is
plot(aso_3m, main = "ASO 3m snow depths (m)")
plot(ER_proj, border = "black", add = TRUE)

### Plot the ASO data as snow/no snow
## Any values zero or below are classified as "no snow" and given a value of 0
## Any values above zero are classified as "snow" and given a value of 1

aso_50m_binary<- reclassify(aso_50m, cbind(-Inf, 0, 0))
aso_50m_binary<- reclassify(aso_50m_binary, cbind(0.1, Inf, 1))
 
aso_3m_binary<- reclassify(aso_3m, cbind(-Inf, 0, 0))
aso_3m_binary<- reclassify(aso_3m_binary, cbind(0.1, Inf, 1))

par(mfrow=c(1,2)) 
plot(aso_50m_binary,  col = rev(cm.colors(40)), main = "ASO 50-m Binary")
plot(ER_proj, border = "black", add = TRUE)
plot(aso_3m_binary, col = rev(cm.colors(40)), main = "ASO 3-m Binary")
plot(ER_proj, border = "black", add = TRUE)


## Per HP, exclude anything <= 20cm due to uncertainties in data below this level
aso_50m_binary<- reclassify(aso_50m, cbind(-Inf, 0.2, 0))
aso_50m_binary<- reclassify(aso_50m_binary, cbind(0.2, Inf, 1))

aso_3m_binary<- reclassify(aso_3m, cbind(-Inf, 0.2, 0))
aso_3m_binary<- reclassify(aso_3m_binary, cbind(0.2, Inf, 1))

### Format and display the categorical aso data
aso_50_bi_fac<- as.factor(aso_50m_binary)
rat_aso<- levels(aso_50_bi_fac)[[1]]
rat_aso[["status"]] <- c("no snow", "snow")
levels(aso_50_bi_fac)<- rat_aso

#dev.new()
#par(mar = c(4,4,4,3.5))
levelplot(aso_50_bi_fac, att = "status", col.regions = rev(cm.colors(40)), scales = list(x = list(cex = 1.5), y = list(cex = 1.5, rot = 90)), xlab = list(label = "Meters", cex = 2), ylab = list(label = "Meters", cex = 2),main = list(label = "ASO Snow Status", cex = 2),  colorkey = list(labels = list(height = 2, cex = 1.9))) + layer(sp.polygons(ER_proj))


### Plot a histogram of the depth values

## Convert rasters into vectors
aso_50m_vec <- as.vector(aso_50m)
aso_3m_vec <- as.vector(aso_3m)

## Plot the original values of each in histograms to visualize their distribution
par(mfrow=c(1,2))
hist(aso_50m_vec, main = "ASO 50m depth values", xlab = "depth (m)", col = "blue")
hist(aso_3m_vec, main = "ASO 3m depth values", xlab = "depth (m)", col = "blue")



### Compare the ASO data to STARFM output from that date
## Bring in the model data and format it to be plotted

starfm<- brick('./starfm_fusion.tif')

## If importing the entire raster brick of STARFM results, specify the layer of interest whose date matches the ASO data
starfm_day1<- brick(starfm[[1]])


### Mask the no data values
starfm_day1<- reclassify(starfm_day1, cbind(-Inf, -11111, NA))
starfm_day1<- starfm_day1/10000


## Need to make the 2 raster projections match

aso_crs = crs(aso_50m)
starfm_day1_proj<- projectRaster(starfm_day1, crs = aso_crs)


### Creating binary snow rasters (with rasters in same projection)

# Any values in our rasters of below 0.4 are given a new value of zero (no snow)
starfm_binary_proj<- reclassify(starfm_day1_proj,cbind(-Inf, 0.4, 0))

# Any values in our rasters of 0.4 or above are given a new value of 1 (snow)
starfm_binary_proj<- reclassify(starfm_binary_proj,cbind(0.4, Inf, 1))

starfm_bi_fac_proj<- as.factor(starfm_binary_proj)
rat_proj<- levels(starfm_bi_fac_proj)[[1]]
rat_proj[["status"]] <- c("no snow", "snow")
levels(starfm_bi_fac_proj)<- rat_proj


levelplot(starfm_bi_fac_proj, col.regions = rev(cm.colors(40)),scales = list(x = list(cex = 1.5), y = list(cex = 1.5, rot = 90)), xlab = list(label = "Meters", cex = 2), ylab = list(label = "Meters", cex = 2),main = list(label = "STARFM Snow Status", cex = 2),  colorkey = list(labels = list(height = 2, cex = 1.9))) + layer(sp.polygons(ER_proj))


## Plot ASO and projected model data together

par(mfrow=c(1,2)) 
plot(aso_50m_binary, col = rev(cm.colors(40)), main = "ASO 50m Data")
plot(ER_proj, border = "black", add = TRUE)
plot(starfm_binary_proj, col = rev(cm.colors(40)), main = "STARFM Data")
plot(ER_proj, border = "black", add = TRUE)

### Find where the data for the two rasters overlaps


## apply the 20cm threshold to the data
aso_cutoff<- aso_50m<- reclassify(aso_50m, cbind(-Inf, 0.2, 0))

## resample the ASO data to 30m pixels
aso_30m <- resample(aso_cutoff, starfm_binary_proj, method = 'bilinear')

## Set all pixels with any data to a value of 1 (regardless of snow status)
aso_30m_data<- reclassify(aso_30m, cbind(-Inf, Inf, 1))

## Using STARFM data where NAs have been assigned and NDSI magnitude corrected, but is not yet projected, set all pixels w/ data to a value of 1
starfm_data<- reclassify(starfm_day1, cbind(-Inf, Inf, 1))

## Now project the STARFM raster
starfm_data_proj<- projectRaster(starfm_data, crs = aso_crs)

## Plot to make sure pixels w/ data were assigned correctly (everything should equal 1)
par(mfrow=c(1,2)) 
plot(aso_30m_data, main = "ASO 30m Data")
plot(ER_proj, border = "black", add = TRUE)
plot(starfm_data_proj, main = "STARFM Data")
plot(ER_proj, border = "black", add = TRUE)

## Use raster math to find overlapping pixels (plots where only both have data)

data_sum <- overlay(aso_30m_data, starfm_data_proj, fun=function(r1,r2){return(r1+r2)})

## Plot only the area where the data overlap
dev.new(height=0.91*nrow(data_sum)/50, width=1.09*ncol(data_sum)/50)
par(mar = c(5,5,5,3.5))
plot(data_sum, col=topo.colors(1), main = "", xlab = "Meters", ylab = "Meters", cex.lab = 2, cex.axis = 1.5, legend = FALSE)
title(main = list("Data Overlap", cex = 2), line = 2.1)
mtext("ASO/STARFM", line = 0.4, cex = 1.7)
plot(ER_proj, border = "black", add = TRUE)

### Consider the entire area of the STARFM raster, calculate the fraction of STARFM data that the ASO data overlaps

## calculate the size of the rasters
size<- nrow(data_sum) * ncol(data_sum)

data_both <- reclassify(data_sum, cbind(1.5, 2.5, 1)) ## set areas with data back to a value of 1

## find the number of pixels with overlapping data in both rasters
data_both_vec <- as.vector(data_both)
both_sum<- sum(data_both_vec, na.rm = TRUE)

## find the fraction of data available for the rasters
total_data<- both_sum/size

print(total_data)

### Consider only the area of ASO data acquisition, calculate the fraction of ASO data that the STARFM data overlaps

aso_30m_vec<- as.vector(aso_30m_data)
aso_sum<- sum(aso_30m_vec, na.rm = TRUE)

total_aso_overlap<- both_sum/aso_sum

print(total_aso_overlap)

### Calculate the accuracy of the model from the ASO data
## Create a binary map for the resampled, 30m data

aso_30m_binary<- reclassify(aso_30m, cbind(-Inf, 0, 0)) ## no snow pixels
aso_30m_binary<- reclassify(aso_30m_binary, cbind(0, Inf, 1)) ## snow pixels

plot(aso_30m_binary, col = rev(cm.colors(40)), main = "ASO 30m Binary")
plot(ER_proj, border = "black", add = TRUE)

## Create a confusion matrix for the model values that overlap the ASO data

aso_30m_binary_vec<- as.vector(aso_30m_binary, mode = 'numeric')
starfm_binary_vec<- as.vector(starfm_binary_proj, mode = 'numeric')

aso_factor<- as.factor(aso_30m_binary_vec)
starfm_factor<- as.factor(starfm_binary_vec)

conf_matrix<- confusionMatrix(starfm_factor,aso_factor, positive = "1")

# Save the results in matrices that can easily be viewed for reference
overall<- as.matrix(conf_matrix, what = "overall")
classes<- as.matrix(conf_matrix, what = "classes")


# Want to know the individual values, so save these values as their own variable
acc_model<- overall[[1]]
spec_model<- classes[[2]]
prec_model<- classes[[5]]
recall_model<- classes[[6]]
f1_model<- classes[[7]]

df<- "ASO Data Comparison"
results_end<- c(df, acc_model, spec_model, recall_model, prec_model, f1_model)
labels<- c("Analysis", "Accuracy", "Specificity", "Recall", "Precision", "F-Score")

all_data<- data.frame(labels, results_end)
