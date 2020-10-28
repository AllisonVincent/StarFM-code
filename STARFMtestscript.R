######## R code for GEE/STARFM ######### 
#### based on code by Faye Peters
#### modified by Megan Gallagher with help from Jake Graham
#### further modification for use with NDSI data by Allison Vincent
######################################################
### Used in conjunction with GEE code, this is for Landsat8 and MODIS Terra Daily 500 m use
### Set up a folder that contains the Landsat, and MODIS files, as well as the starfm.exe and starfmconfig.txt
### create a separate folder called output in that folder

# install.packages("rdgal", "fields", "sp", "zoo", "ggplot2", "raster", "spam", "maps")

# modified with help from Jake on 6/10/2019


library(raster)
library(rgdal)
library(zoo)
library(ggplot2)
library(fields)

## Load images and dates (be careful of matching dates and layers, the first layer may repeat)
## Read imagery dates and matching data:
## image_dates<- read.csv("./2016_Dates_BOP_East_run1.csv", stringsAsFactors=F) 
## dates are not used in this code but are useful for other programs like bfastspatial, to check date overlap, and timing

modis<-brick("./2016_mod_East.tif")
landsat<-brick("./2016_landsat_East.tif")

##percentage of cover that is used to find "good" landsat scenes for interpolation 
## is the inverse of the actual cloud cover, 99 percent means 1 percent of the pixels is an actual value
perc<- 99

##################################################Code########################################################
#### STARFM code and function

# STARFM Function
 #  INPUTS:
 #     modis = raster brick with MODIS data for all dates
 #     landsat = raster brick with all Landsat data for time period, also nodata fill layers for dates w/out Landsat
 #     perc = percentage of cover used to find "good" landsat scenes for interpolation
 #     handleNA = logical variable indicating if the zeros in the data have been replaced with NA yet, default to TRUE for the first run

 #  OUTPUTS: all outputs (except landsat_sim) updated with each loop iteration
 #     "./modis_t1.envi" data for modis at the date corresponding with the "before" landsat date 
 #     "./modis_t2.envi" data for modis, regardless of landsat status 
 #     "./modis_t3.envi" data for modis at the date corresponsing with the "after" landsat date 
 #     "./landsat_t1.envi" data for landsat at the "before" date
 #     "./landsat_t3.envi" data for landsat at the "after" dates
 #     "./landsat_t2_sim.envi" simulated data 
 #     landsat_sim simulated data with NA values


# reads the number of rows and columns in the input rasters and pastes their values into the config file
starfm<-function(modis,landsat,perc, handleNA = TRUE){
  config <- readLines("./StarFM_config.txt")
  config <- gsub("(.*NROWS = ).*$", paste0("\\1", nrow(landsat)),config )
  config <- gsub("(.*NCOLS = ).*$", paste0("\\1", ncol(landsat)),config )
  cat(config, file="./StarFM_config.txt", sep="\n")
  
  

  if(handleNA){
     ## Fix any missing data that enters as zeroes to our NA value. We have
     ## to do this by layer as operating on the entire stack may run into
     ## memory issues on larger subsets:
     for (i in 1:nlayers(modis)) {
       ## Use the raster Which() function for speed:
       masked <- Which(modis[[i]] == 0, cells=TRUE)
       modis[[ i ]][ masked ] <- -32768
       masked <- Which(landsat[[i]] == 0, cells=TRUE)
       landsat[[ i ]][ masked ] <- -32768
     }
    writeRaster(modis, "./2016_mod_NA_Handled.tif")
    writeRaster(landsat,"./2016_landsat_NA_Handled.tif")
  } else {
    masked <- brick("./2016_mod_NA_Handled.tif")
    landsat <- brick("./2016_landsat_NA_Handled.tif")
  }
  print("here 2")
  flush.console()
  
  ## Automatically choose "good" landsat layers, at the moment this includes all actual landsat images
  ## if you want to create a threshhold for masking change the total percent of the pixel.
  ## i.e. if percent = 30, more than 70 percent of the image must be actual pixel values
  
  
  area<-landsat@nrows*landsat@ncols  # find total area of landsat grid
  perc_Area<-(perc/100)*area  # find the value of (example 99%) total grid area
  
  # Preallocating variables
  test2<-rep(NA,nlayers(modis)) # create a new logical vector with number of NA values as raster brick layers
  filternew<-rep(NA,nlayers(modis))  # same here, only new variable name
  
  
  # Go through landsat layers (starting at the 2nd layer). Assuming that the first layer is a valid landsat layer  
  for (i in 2:nlayers(modis)){  
    test2[i] <-sum(landsat[[i]][])  # sum the value of the data found in each layer 
    filternew[i] <-(test2[i]>{-32768*perc_Area}) == 1  # if the sum from above is greater than -32768*perc_Area, then set the value of that layer equal to TRUE
  }
  
  
  filt3 = which(filternew==1) # find the positions of the layers above equal TRUE
  filt3 <-append(filt3,1,0) # add a 1 to the beginning of the filt3 vector
  # end up with the filt3 variable being the layer id numbers of all legitimate landsat data (every Landsat date included here, the layers in between each 16 days removed)
  
  ## Test a single landsat image from list
  
  good_layer <-filt3[3] # select the 3rd landsat image (could be any image)
  plot(landsat[[good_layer]]) # plot it to look at it
  
  ## If the above all works, then we run the following to loop over the
  ## MODIS time steps, filling in Landsat output as we go:
  
  landsat_sim <- stack(modis) # duplicate the modis stack as a new variable
  landsat_sim[] <- NA # convert all values in the stack to NA
  
  ## Iterate and run StarFM for each MODIS date, choosing the
  ## nearest pair of good MODIS/Landsat dates, one before and
  ## one after the date being simulated where possible: 
  print("HERE!")
  flush.console()
  
  good_landsat <- c(filt3) # create a new vector with the positions of the valid landsat layers from above
  if(!length(good_landsat)){ # if good_landsat has no values (i.e. is NULL0
    print("WARNING!!!! No good landsat dates... exiting program...")
    flush.console()
    return(-1)
  } else {
    pb <- pbCreate((nlayers(landsat_sim)), "window", style=3,label='Time Step Progress') # create a progress bar
    for (i in 1:nlayers(landsat_sim)) { # for each layer in landsat sim rasterbrick
      
      ## jakes/megan's work
      if (i %in% c(filt3)){  #if the layer in the landsat sim rasterbrick is also in the vector filt3.. i.e., you have a 'good' landsat image on this day
        ls_t1<- i # set variable equal to layer id.. there is a "good" landsat image on this day, so no need to pull post/prev images
        ls_t3 <- i # set variable equal to layer id.. there is a "good" landsat image on this day, so no need to pull post/prev images
      }
      else {  # if the layer (i.e., day) does not have a "good" landsat image then do...
        foo <- good_landsat - i # subtract i (~DOY) from the "good" landsat dates. This produces a measure of time differences
        if(!length(which(foo==0))){ # safeguard in case there are NO "good" landsat images. Is there a value in foo equal to zero?
          ls_t1 <- good_landsat[which(foo==max(foo[foo < 0]))] # select the closest date with a "good" image BEFORE this date
          ls_t3  <- good_landsat[which(foo==min(foo[foo > 0]))] # select the closest date with a "good" image AFTER this date
        }
      }
      
      
      ##end
      
      m_t1 <- ls_t1  # set "good" landsat date to a variable that can be used for modis "before" pair
      m_t3 <- ls_t3  # set "good" landsat date to a variable that can be used for modis "after" pair

      modis_t1 <- modis[[m_t1]]  # raster of modis at date with corresponding "before" landsat dates
      modis_t2 <- modis[[i]]  # raster of modis at dates at that step, regardless of landsat  
      modis_t3 <- modis[[m_t3]] # raster of modis at date with corresponding "after" landsat dates
      landsat_t1 <- landsat[[ls_t1]]  # raster of good landsat at "before" date
      landsat_t3 <- landsat[[ls_t3]]  # raster of good landsat at "after" date
      
      ## write rasters for StarFM to work on (re-written for every step) ... can be seen in config file... "StarFM_config.txt"
      writeRaster(modis_t1, filename="./modis_t1.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
      writeRaster(modis_t2, filename="./modis_t2.envi",bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
      writeRaster(modis_t3, filename="./modis_t3.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
      writeRaster(landsat_t1, filename="./landsat_t1.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
      writeRaster(landsat_t3, filename="./landsat_t3.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
      
      # run the actual STARFM model
      system2(command="./StarFM.exe",args="StarFM_config.txt", wait=TRUE) # Produces "/landsat_t2_sim.envi" 
      
      landsat_t2_sim <- raster("./landsat_t2_sim.envi")
      
      ## Set any -32768 to NA values before writing:
      landsat_t2_sim[ landsat_t2_sim == -32768 ] <- NA
      landsat_sim[[i]] <- landsat_t2_sim[]
      
      ## In our filled data set, set any missing Landsat pixels to those
      ## simulated via StarFM:
      
      
    pbStep(pb, step=NULL, label='Processed Layer') } # progress bar related..
    pbClose(pb, timer=T) # progress bar.. 
    return(landsat_sim)
  }
}

#### here we actually run the function

#first time running
landsat_sim<-starfm(modis,landsat,perc)

#all other times, uncomment if input rasters with NA values have already been created
# landsat_sim<-starfm(modis,landsat,perc, handleNA = F)


#### Possible raster outputs for saving
## writeRaster(landsat_sim, filename="./output/2016_East_run1_fusion.envi",
##            bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)

writeRaster(landsat_sim, filename="./output/2016_East_fusion.tif", bandorder='BSQ', 
            datatype='INT2S',format='GTiff', overwrite=TRUE)





