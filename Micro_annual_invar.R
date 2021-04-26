#### run parts of Fridley microclimate model that are the same every year

#### setup ####
  # packages
    library(raster)
    library(rgdal)
    library(lubridate)
    library(parallel)


  # paths
    fridley_GIS <- "E:/Smokies_Veg/phys_GIS/gsmnp_ascii/"
    climate_path <- "E:/Smokies_Veg/climate_data/"  
  
  # import GIS data
    tci <- log(raster(paste(fridley_GIS,"tci_cor.asc",sep="")))
    totrad <- raster(paste(fridley_GIS,"totrad.txt",sep=""))
    elev <- raster(paste(fridley_GIS,"elev.txt",sep=""))
    strdist <- raster(paste(fridley_GIS,"logsd.txt",sep=""))
    
    radnames <- list.files(paste(fridley_GIS,"rad/",sep=""),full.names=T)
    
    radstack <- stack(radnames,quick=T)

    
    
    crs(radstack) <- crs(elev) <- crs(strdist) <- crs(totrad) <- crs(tci) <- CRS("+proj=utm +zone=17 +datum=NAD27")
    
    tci <- crop(tci,radstack)
    totrad <- crop(totrad,radstack)
    elev <- crop(elev,radstack)  
    strdist <- crop(strdist,radstack)
    
    names(tci) <- "tci"
    names(totrad) <- "totrad"
    names(elev) <- "elev"
    names(strdist) <- "strdist"
      
  # import model coefficients
    maxtmod <- read.csv(paste(climate_path,"micro_mod/maxcoef_out.csv",sep=""),row.names=1)
    mintmod <- read.csv(paste(climate_path,"micro_mod/mincoef_out.csv",sep=""),row.names=1)

## run calculation in parallel and save files
    CalcTMP <- function(doy){
      
      rad <- radstack[[doy]]
      jday <- doy
      cosday <- cos(0.0172 * jday)
      sinday <- sin(0.0172 * jday)
      
      tmp_min <- mintmod["(Intercept)",1] + 
        (mintmod["ELEV",1] * elev) +
        (mintmod["LOG.STRDST",1] * strdist) +
        (mintmod["I(log(TCI))",1] * tci) +
        (mintmod["RAD",1] * rad) +
        (mintmod["cos(0.0172 * JDATE)",1] * cosday)+ 
        (mintmod["sin(0.0172 * JDATE)",1] * sinday) +
        (mintmod["cos(0.0172 * JDATE):ELEV",1] * cosday * elev) + 
        (mintmod["sin(0.0172 * JDATE):ELEV",1] * sinday * elev) +
        (mintmod["cos(0.0172 * JDATE):LOG.STRDST",1] * cosday * strdist)  + 
        (mintmod["sin(0.0172 * JDATE):LOG.STRDST",1] * sinday * strdist) +
        (mintmod["cos(0.0172 * JDATE):I(log(TCI))",1] * cosday * tci) + 
        (mintmod["sin(0.0172 * JDATE):I(log(TCI))",1] * sinday * tci) +
        (mintmod["RAD:ELEV",1] * rad * elev) +
        (mintmod["RAD:LOG.STRDST",1] * rad * strdist)
      
      writeRaster(tmp_min,
                  paste(climate_path,"tmpMin/d_",sprintf("%03d",doy),sep=""),
                  format = "GTiff")
      
      tmp_max <- maxtmod["(Intercept)",1]  + 
        (maxtmod["TOTRAD",1] * totrad) +
        (maxtmod["ELEV",1] * elev) +
        (maxtmod["LOG.STRDST",1] * strdist) +
        (maxtmod["I(log(TCI))",1] * tci) +
        (maxtmod["RAD",1] * rad) +
        (maxtmod["cos(0.0172 * JDATE)",1] * cosday) + 
        (maxtmod["sin(0.0172 * JDATE)",1] * sinday) +
        (maxtmod["cos(0.0172 * JDATE):TOTRAD",1] * cosday * totrad) +  
        (maxtmod["sin(0.0172 * JDATE):TOTRAD",1] * sinday * totrad) +
        (maxtmod["cos(0.0172 * JDATE):ELEV",1] * cosday * elev)   + 
        (maxtmod["sin(0.0172 * JDATE):ELEV",1] * sinday * elev) +
        (maxtmod["cos(0.0172 * JDATE):LOG.STRDST",1] * cosday * strdist) + 
        (maxtmod["sin(0.0172 * JDATE):LOG.STRDST",1] * sinday * strdist) +
        (maxtmod["cos(0.0172 * JDATE):I(log(TCI))",1] * cosday * tci)    + 
        (maxtmod["sin(0.0172 * JDATE):I(log(TCI))",1] * sinday * tci) +
        (maxtmod["RAD:ELEV",1] * rad * elev) +
        (maxtmod["RAD:I(log(TCI))",1] * rad * tci)
      
      writeRaster(tmp_max,
                  paste(climate_path,"tmpMax/d_",sprintf("%03d",doy),sep=""),
                  format="GTiff")
      
      
      
    }
    
    
    clus <- makePSOCKcluster(15)
    clusterExport(cl=clus,varlist=ls())
    clusterEvalQ(cl=clus,library(raster))
    clusterApply(cl=clus,x=1:365,fun=CalcTMP)
    stopCluster(clus)
    
    
    # avoid writng GIS variables back to global env
    allvars <- ls()
    
    rmvars <- allvars[!allvars %in% mainvars]
    
    rm(list=rmvars)
    