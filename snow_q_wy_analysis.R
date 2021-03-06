### This script is for plotting annual stream flow and SCA (in percent of the watershed that is snow covered) on the same graph

### Stream flow data is from the USGS (water.usgs.gov)


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

## Rasters of the total number of snow covered days per pixel by water year 
snow_sum_2008<- raster("./WY2008_snow_sum.tif")
snow_sum_2010<- raster("./WY2010_snow_sum.tif")
snow_sum_2012<- raster("./WY2012_snow_sum.tif")

## Project the watershed shapefile to match the projection of the raster data
data_proj<- crs(snow_sum_2008)
ER_proj<- spTransform(ER, data_proj)
proj4string(ER_proj)

## Rasters of the total number of days each pixel had data (regardless of snow status) by water year
total_model_2008<- raster('./WY2008_model_sum.tif')
total_model_2010<- raster('./WY2010_model_sum.tif')
total_model_2012<- raster('./WY2012_model_sum.tif')


### Using the raster data above, find the percent of snow-covered days by pixel for each water year
## Start with WY 2008
corrected_2008<- (snow_sum_2008/total_model_2008) * 100

## Turn the above into a vector and find the mean and median of the data
corrected_2008_vec<- as.vector(corrected_2008)
mean_2008<- mean(corrected_2008_vec)
median_2008<- median(corrected_2008_vec)
corrected_2008_df<- as.data.frame(corrected_2008_vec)


## Plot the per-pixel percent snow-cover as a histogram
hist_2008<- ggplot(corrected_2008_df, aes(x = corrected_2008_vec)) + 
  geom_histogram(bins = 30, 
                 col = "black",
                 fill = "blue") +
  ylim(c(0, 500000)) +  ## this parameter allows for setting a common y-axis limit for all plots
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15)) +
  labs(title = "2008 Percent of Year Pixels are Snow Covered", x = "Percent snow covered", y = "Frequency")

hist_2008


## Option to plot the mean and median values as vertical lines on the histogram
#hist_2008 +
#  geom_vline(aes(xintercept = mean(corrected_2008_vec)), col = "red", size = 1) +
#  geom_text(aes(label = round(mean_2008, 1), y = 0, x = mean_2008), vjust = -1, col = "red", size = 4) +
#  geom_vline(aes(xintercept = median(corrected_2008_vec)), col = "purple", size = 1) 
#  geom_text(aes(label = round(median_2008, 1), y = 0, x = median_2008), vjust = -1, col = "purple", size = 4)


### Repeat the above, but for WY 2010
corrected_2010<- (snow_sum_2010/total_model_2010) * 100
corrected_2010_vec<- as.vector(corrected_2010)
corrected_2010_df<- as.data.frame(corrected_2010_vec)


hist_2010<- ggplot(corrected_2010_df, aes(x = corrected_2010_vec)) +
  geom_histogram(bins = 30,
                 col = "black",
                 fill = "blue") +
  ylim(c(0, 500000)) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15)) +
  labs(title = "2010 Percent of Year Pixels are Snow Covered", x = "Percent snow covered", y = "Frequency")

hist_2010

# hist_2010 +
#   geom_vline(aes(xintercept = mean(corrected_2010_vec)), col = "red", size = 1) +
#   geom_vline(aes(xintercept = median(corrected_2010_vec)), col = "purple", size = 1)


### Repeat the above, but for WY 2012
corrected_2012<- (snow_sum_2012/total_model_2012) *100
corrected_2012_vec<- as.vector(corrected_2012)
corrected_2012_df<- as.data.frame(corrected_2012_vec)

hist_2012<- ggplot(corrected_2012_df, aes(x = corrected_2012_vec)) +
  geom_histogram(bins = 30,
                 col = "black",
                 fill = "blue") +
  ylim(c(0, 500000)) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15)) +
  labs(title = "2012 Percent of Year Pixels are Snow Covered", x = "Percent snow covered", y = "Frequency")

hist_2012

# hist_2012 +
#   geom_vline(aes(xintercept = mean(corrected_2012_vec)), col = "red", size = 1) +
#   geom_vline(aes(xintercept = median(corrected_2012_vec)), col = "purple", size = 1)


### Load the table that has the total number of snow-covered pixels by date and percentage of the watershed that is snow covered

wy2008_total<- read.csv("./WY2008/WY2008_snow_sum_by_day.csv")
wy2010_total<- read.csv("./WY2010/WY2010_snow_sum_by_day.csv")
wy2012_total<- read.csv("./WY2012/WY2012_snow_sum_by_day.csv")

## format the dates in the above tables
wy2008_total$Date<- as.Date(wy2008_total$Date, format = "%m/%d/%Y")
wy2010_total$Date<- as.Date(wy2010_total$Date, format = "%m/%d/%Y")
wy2012_total$Date<- as.Date(wy2012_total$Date, format = "%m/%d/%Y")


### Plot just the percent of the watershed that is snow covered for each water year
ggplot(data = wy2008_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2008")


ggplot(data = wy2010_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2010")


ggplot(data = wy2012_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2012")

### Try plotting all of the above on one graph just to see what it looks like

## Load a file that has the percent snow-covered data for all water years in one
perc_all_wy<- read.csv(file = './inter_annual/snow_perc_all_years.csv')


DOWY<- perc_all_wy$DOWY # DOWY = 'day of water year'
wy2008<- perc_all_wy$X2008_perc  
wy2010<- perc_all_wy$X2010_perc
wy2012<- perc_all_wy$X2012_perc


ggplot(perc_all_wy) +
  geom_line(aes(x = DOWY, y = wy2008, color = "blue")) +
  geom_line(aes(x = DOWY, y = wy2010, color = "black")) +
  geom_line(aes(x = DOWY, y = wy2012, color = "red")) +
  scale_color_identity(name = "Water Year",
                       breaks = c("blue", "black", "red"),
                       labels = c("WY 2008", "WY 2010", "WY 2012"),
                       guide = "legend") +
  ggtitle("Perc snow covered pixels by water year") +
  xlab("Day of Water Year") +
  ylab("Percent")


##### Use the function below to find the moving average for each water year

# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {
  
  
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

#### Plot the daily percent snow values as bars and the moving average as a line on top

wy_2008_avg<- movingAverage(wy2008, 10, TRUE)   ## A 10-day centered moving average seems to be the minimum for a smoother curve
perc_all_wy$mvg_avg_2008<- wy_2008_avg

ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2008, fill = "sky blue"), color = "sky blue") + 
  geom_line(aes(x = DOWY, y = wy_2008_avg, color = "firebrick"), size = 1.0) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("Perc pixels snow covered"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick"),
                       labels = c("Moving Avg"),
                       guide = "legend") +
  ggtitle("Perc snow covered pixels for 2008 water year") +
  xlab("Day of Water Year") +
  ylab("Percent")


wy_2010_avg<- movingAverage(wy2010, 10, TRUE)
perc_all_wy$mvg_avg_2010<- wy_2010_avg

ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2010, fill = "sky blue"), color = "sky blue") + 
  geom_line(aes(x = DOWY, y = wy_2010_avg, color = "firebrick"), size = 1.0) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("Perc pixels snow covered"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick"),
                       labels = c("Moving Avg"),
                       guide = "legend") +
  ggtitle("Perc snow covered pixels for 2010 water year") +
  xlab("Day of Water Year") +
  ylab("Percent")


wy_2012_avg<- movingAverage(wy2012, 10, TRUE)
perc_all_wy$mvg_avg_2012<- wy_2012_avg

ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2012, fill = "sky blue"), color = "sky blue") + 
  geom_line(aes(x = DOWY, y = wy_2012_avg, color = "firebrick"), size = 1.0) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("Perc pixels snow covered"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick"),
                       labels = c("Moving Avg"),
                       guide = "legend") +
  ggtitle("Perc snow covered pixels for 2012 water year") +
  xlab("Day of Water Year") +
  ylab("Percent")



#### Plot all the moving average lines above on the same graph to compare them


ggplot(perc_all_wy) +
  geom_line(aes(x = DOWY, y = wy_2008_avg, color = "blue"), size = 1.0) +
  geom_line(aes(x = DOWY, y = wy_2010_avg, color = "black"), size = 1.0) +
  geom_line(aes(x = DOWY, y = wy_2012_avg, color = "red"), size = 1.0) +
  scale_color_identity(name = "Water Year",
                       breaks = c("blue", "black", "red"),
                       labels = c("WY 2008", "WY 2010", "WY 2012"),
                       guide = "legend") +
  ggtitle("Moving avg snow covered pixels by water year") +
  xlab("Day of Water Year") +
  ylab("Percent")

