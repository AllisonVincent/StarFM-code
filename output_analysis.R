### This script is for viewing and acquiring basic information about individual layers of Landsat, MODIS, and STARFM data


library(raster)

landsat<- brick('./landsat.tif')
modis<- brick('./mod_.tif')

data<- brick('./starfm_East_fusion.tif') ## starfm data



## To find the fraction of data available, first set pixels with data to 1, then convert all non-data values for each raster layer to NA

## First, do this for STARFM results

for (i in 1:nlayers(data)) {
  # Use the raster Which() function for speed:
  masked <- Which(data[[i]] > -11111 , cells=TRUE)
  data[[i]][masked] <- 1
}


for (i in 1:nlayers(data)) {
  # Use the raster Which() function for speed:
  masked <- Which(data[[i]] == -11111, cells=TRUE)
  data[[i]][masked] <- NA
}


## Compute the data fraction for each layer

## total number of pixels in raster
area<- ncol(data) * nrow(data)

n <- nlayers(data)

data_frac<- rep(NA, n) ## create an empty raster with the required number of layers to fill in via the loop below:

for (i in 1:nlayers(data)) { 
  layer<- as.vector(data[[i]], mode = 'numeric')
  good_data<- sum(layer, na.rm = TRUE) ## find the sum, or the number of pixels with data, for each layer
  frac<- good_data/area ## calculate the fraction of pixels with data for each layer
  data_frac[[i]]<- frac
}

starfm_df<- data.frame("STARFM" = data_frac) ## convert above results into a data frame to write to a .csv table file

#write.table(starfm_df, "./starfm_data.csv", row.names = FALSE)


################## Repeat the above, but for Landsat


for (i in 1:nlayers(landsat)) {
  # Use the raster Which() function for speed:
  masked <- Which(landsat[[i]] == 0 , cells=TRUE)
  landsat[[i]][masked] <- NA
}


for (i in 1:nlayers(landsat)) {
  # Use the raster Which() function for speed:
  masked <- Which(landsat[[i]] != 0 , cells=TRUE)
  landsat[[i]][masked] <- 1
}


## Compute the data fraction for each landsat layer

area<- ncol(landsat) * nrow(landsat)

n <- nlayers(landsat)

land_frac<- rep(NA, n)

for (i in 1:nlayers(landsat)) { 
  layer<- as.vector(landsat[[i]], mode = 'numeric')
  good_data<- sum(layer, na.rm = TRUE)
  frac<- good_data/area
  land_frac[[i]]<- frac
}


################## Now do the same for modis data


for (i in 1:nlayers(modis)) {
  # Use the raster Which() function for speed:
  masked <- Which(modis[[i]] != 0 , cells=TRUE)  ## no data values for MODIS area already set to NA
  modis[[i]][masked] <- 1
}

## Compute the data fraction for each layer

area<- ncol(modis) * nrow(modis)

n <- nlayers(modis)

modis_frac<- rep(NA, n)

for (i in 1:nlayers(modis)) { 
  layer<- as.vector(modis[[i]], mode = 'numeric')
  good_data<- sum(layer, na.rm = TRUE)
  frac<- good_data/area
  modis_frac[[i]]<- frac
}


## Put all the above data into a single dataframe 
df<- data.frame("Landsat" = land_frac, "Modis" = modis_frac, "STARFM" = data_frac)
