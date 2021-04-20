setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Analysis_Tests/Stream_data')

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

### Read in the csv file

df_2008<- read.csv(file = 'wy_2008.csv')
df_2010<- read.csv(file = 'wy_2010.csv')
df_2012<- read.csv(file = 'wy_2012.csv')


### Convert dates from characters to dates

df_2008$Date<- as.Date(df_2008$Date, format = "%m/%d/%Y")
df_2010$Date<- as.Date(df_2010$Date, format = "%m/%d/%Y")
df_2012$Date<- as.Date(df_2012$Date, format = "%m/%d/%Y")


### Map discharge from individual water years

ggplot(data = df_2008, mapping = aes(x = Date, y = cfs)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Discharge (cfs)",
    title = "Water Year 2008"
  )


ggplot(data = df_2010, mapping = aes(x = Date, y = cfs)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Discharge (cfs)",
    title = "Water Year 2010"
  )


ggplot(data = df_2012, mapping = aes(x = Date, y = cfs)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Discharge (cfs)",
    title = "Water Year 2012"
  )

###### Plot all discharge data together

df_all<- read.csv(file = 'all_wys.csv')


DOWY<- df_all$DOWY
wy2008_q<- df_all$X2008
wy2010_q<- df_all$X2010
wy2012_q<- df_all$X2012


ggplot(df_all) +
  geom_line(aes(x = DOWY, y = wy2008_q, color = "blue")) +
  geom_line(aes(x = DOWY, y = wy2010_q, color = "black")) +
  geom_line(aes(x = DOWY, y = wy2012_q, color = "red")) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 15)) +
  scale_color_identity(name = "Water Year",
                       breaks = c("blue", "black", "red"),
                       labels = c("WY 2008", "WY 2010", "WY 2012"),
                       guide = "legend") +
  ggtitle("Discharge at Watershed Outlet by Water Year") +
  xlab("Day of Water Year") +
  ylab("Discharge (cfs)")


##### Create cumulative discharge plots for individual water years

df_2008[, "cum_Q"] <- cumsum(df_2008$cfs)
df_2010[, "cum_Q"] <- cumsum(df_2010$cfs)
df_2012[, "cum_Q"] <- cumsum(df_2012$cfs)


ggplot(data = df_2008, mapping = aes(x = Date, y = cum_Q)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Cumulative Discharge (cfs)",
    title = "Water Year 2008"
  )


ggplot(data = df_2010, mapping = aes(x = Date, y = cum_Q)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Cumulative Discharge (cfs)",
    title = "Water Year 2010"
  )


ggplot(data = df_2012, mapping = aes(x = Date, y = cum_Q)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Cumulative Discharge (cfs)",
    title = "Water Year 2012"
  )


####### Create cumulative sum discharge plot with all water years

## Two of our WYs are leap years, so they contain one extra value than WY 2010. To do a cumulative sum, we need to replace the NA for this day with a value of 0

temp_wy2010<- df_all$X2010
temp_wy2010[is.na(temp_wy2010)] = 0


df_all[, "cum_Q_2008"] <- cumsum(df_all$X2008)
df_all[, "cum_Q_2010"] <- cumsum(temp_wy2010)
df_all[, "cum_Q_2012"] <- cumsum(df_all$X2012)

cumsum_2008<- df_all$cum_Q_2008
cumsum_2010<- df_all$cum_Q_2010
cumsum_2012<- df_all$cum_Q_2012


ggplot(df_all) +
  geom_line(aes(x = DOWY, y = cumsum_2008, color = "blue")) +
  geom_line(aes(x = DOWY, y = cumsum_2010, color = "black")) +
  geom_line(aes(x = DOWY, y = cumsum_2012, color = "red")) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 13)) +
  scale_color_identity(name = "Water Year",
                       breaks = c("blue", "black", "red"),
                       labels = c("WY 2008", "WY 2010", "WY 2012"),
                       guide = "legend") +
  ggtitle("Discharge Cumulative Sums at Watershed Outlet by Water Year") +
  xlab("Day of Water Year") +
  ylab("Discharge (cfs)")

### Find the center of mass, or when 50% of discharge has passed (doing this manually because I'm not sure how else to do it)

com_2008<- cumsum_2008[366]/2  ## according to the cumsum values, 50% of discharge passes through on day 253
com_2010<- cumsum_2010[366]/2  ## according to the cumsum values, 50% of discharge passes through on day 246 (subtract 1 for the leap year)
com_2012<- cumsum_2012[366]/2 ## according to the cumsum values, 50% of discharge passes through on day 221


########### Combine the percent pixel snow covered data with the hydrographs for their respective years

### Load the table that has the total number of snow-covered pixels by date and percentage of the watershed that is snow covered

setwd('C:/Users/Allison and Brian/Documents/Research/STARFM/STARFMtest/Analysis_Tests')

## load the snow percent by day data

wy2008_total<- read.csv("./WY2008/WY2008_snow_sum_by_day.csv")
wy2010_total<- read.csv("./WY2010/WY2010_snow_sum_by_day.csv")
wy2012_total<- read.csv("./WY2012/WY2012_snow_sum_by_day.csv")

wy2008_total$Date<- as.Date(wy2008_total$Date, format = "%m/%d/%Y")
wy2010_total$Date<- as.Date(wy2010_total$Date, format = "%m/%d/%Y")
wy2012_total$Date<- as.Date(wy2012_total$Date, format = "%m/%d/%Y")


ggplot(data = wy2008_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2008"
  )


ggplot(data = wy2010_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2010"
  )


ggplot(data = wy2012_total, mapping = aes(x = Date, y = perc.of.watershed)) +
  geom_line() +
  labs(
    x = "Date",
    y = "Perc of snow covered pixels",
    title = "Water Year 2012"
  )



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

perc_all_wy<- read.csv(file = './inter_annual/snow_perc_all_years.csv')


DOWY<- perc_all_wy$DOWY
wy2008_perc<- perc_all_wy$X2008_perc
wy2010_perc<- perc_all_wy$X2010_perc
wy2012_perc<- perc_all_wy$X2012_perc


### create plot for WY 2008

wy_2008_avg<- movingAverage(wy2008_perc, 10, TRUE)   ## A 10-day centered moving average seems to be the minimum for a smoother curve
perc_all_wy$mvg_avg_2008<- wy_2008_avg

### find the latest day of the water year when percent snow cover is 50% or greater

which(wy_2008_avg >= 50) ## gives all position where the value is above 50%. Take the last one in the array


ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2008_perc, fill = "sky blue"), color = "sky blue") +  ## snow cover percent
  geom_line(aes(x = DOWY, y = wy_2008_avg, color = "firebrick"), size = 1.0) + ## moving avg line
  geom_line(aes(x = DOWY, y = wy2008_q/10, color = "black"), size = 1.0) + ## stream discharge 
  geom_vline(aes(xintercept = 253, color = "purple"), size = 1) +  ## the line that indicates the center of mass
  geom_vline(aes(xintercept = 208, color = "dark green"), size = 1) + ## the line that indicates the last day when the mvg avg of snow cover percent is 59% or greater
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 13)) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("SCA(%)"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick", "black", "purple", "dark green"),
                       labels = c("SCA Moving Avg", "Stream Discharge", "Q center of mass", "Last day of 50% SCA"),
                       guide = "legend") +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2008") +
  xlab("Day of Water Year") 



### create plot for WY 2010

temp_wy2010<- wy2010_perc
temp_wy2010[is.na(temp_wy2010)] = 0

wy_2010_avg<- movingAverage(temp_wy2010, 10, TRUE)
perc_all_wy$mvg_avg_2010<- wy_2010_avg

### find the latest day of the water year when percent snow cover is 50% or greater

which(wy_2010_avg >= 50) ## gives all position where the value is above 50%. Take the last one in the array



ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2010_perc, fill = "sky blue"), color = "sky blue") +
  geom_line(aes(x = DOWY, y = wy_2010_avg, color = "firebrick"), size = 1.0) +
  geom_line(aes(x = DOWY, y = wy2010_q/10, color = "black"), size = 1.0) +
  geom_vline(aes(xintercept = 246, color = "purple"), size = 1) +
  geom_vline(aes(xintercept = 120, color = "dark green"), size = 1) +
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 13)) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("SCA(%)"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick", "black", "purple", "dark green"),
                       labels = c("SCA Moving Avg", "Stream Discharge",  "Q center of mass", "Last day of 50% SCA"),
                       guide = "legend") +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2010") +
  xlab("Day of Water Year") 


### create plot for WY 2012

wy_2012_avg<- movingAverage(wy2012_perc, 30, TRUE)
perc_all_wy$mvg_avg_2012<- wy_2012_avg

which(wy_2012_avg >= 50) ## gives all position where the value is above 50%. Take the last one in the array

 
ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2012_perc, fill = "sky blue"), color = "sky blue") +
  geom_line(aes(x = DOWY, y = wy_2012_avg, color = "firebrick"), size = 1.0) +
  geom_line(aes(x = DOWY, y = wy2012_q/10, color = "black"), size = 1.0) +
  geom_vline(aes(xintercept = 221, color = "purple"), size = 1) +
  geom_vline(aes(xintercept = 165, color = "dark green"), size = 1) +
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 13)) +
  scale_fill_identity(name = NULL,
                      breaks = c("sky blue"),
                      labels = c("SCA(%)"),
                      guide = "legend") +
  scale_color_identity(name = NULL,
                       breaks = c("firebrick", "black", "purple", "dark green"),
                       labels = c("SCA Moving Avg", "Stream Discharge", "Q center of mass", "Last day of 50% SCA"),
                       guide = "legend") +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2012") +
  xlab("Day of Water Year") 



#### Plot the above, but without the legend

ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2008_perc), color = "sky blue") +  ## snow cover percent
  geom_line(aes(x = DOWY, y = wy_2008_avg), color = "firebrick", size = 1.0) + ## moving avg line
  geom_line(aes(x = DOWY, y = wy2008_q/10), color = "black", size = 1.0) + ## stream discharge 
  geom_vline(aes(xintercept = 253), color = "purple", size = 1) +  ## the line that indicates the center of mass
  geom_vline(aes(xintercept = 208), color = "dark green", size = 1) + ## the line that indicates the last day when the mvg avg of snow cover percent is 59% or greater
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13)) +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2008") +
  xlab("Day of Water Year") 


ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2010_perc), color = "sky blue") +
  geom_line(aes(x = DOWY, y = wy_2010_avg), color = "firebrick", size = 1.0) +
  geom_line(aes(x = DOWY, y = wy2010_q/10), color = "black", size = 1.0) +
  geom_vline(aes(xintercept = 246), color = "purple", size = 1) +
  geom_vline(aes(xintercept = 120), color = "dark green", size = 1) +
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13)) +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2010") +
  xlab("Day of Water Year") 



ggplot(perc_all_wy) +
  geom_col(aes(x = DOWY, y = wy2012_perc), color = "sky blue") +
  geom_line(aes(x = DOWY, y = wy_2012_avg), color = "firebrick", size = 1.0) +
  geom_line(aes(x = DOWY, y = wy2012_q/10), color = "black", size = 1.0) +
  geom_vline(aes(xintercept = 221), color = "purple", size = 1) +
  geom_vline(aes(xintercept = 165), color = "dark green", size = 1) +
  scale_y_continuous(
    name = "Percent",
    sec.axis = sec_axis(~.*10, name = "Discharge (cfs)")
  ) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 13)) +
  ggtitle("Percent Snow Cover and Stream Discharge for WY 2012") +
  xlab("Day of Water Year") 
