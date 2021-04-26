#### run Fridley microclimate model; adjust climate predictions for change in lapse

#### setup ####
  # packages
    library(raster)
    library(rgdal)
    library(lubridate)
    library(parallel)

  # paths
    #fridley_GIS <- "E:/Smokies_Veg/phys_GIS/gsmnp_ascii/"
    #climate_path <- "E:/Smokies_Veg/climate_data/"  
  
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
        
    #stack <- stack(tci,totrad,elev,strdist,radstack)
    
  # import lapse rates
    
    maxt_lapse <- read.csv(paste(climate_path,"maxtemp_synoptic_1900_2011.csv",sep=""))
    mint_lapse <- read.csv(paste(climate_path,"mintemp_synoptic_1900_2011.csv",sep=""))
    
    maxt_lapse$date <- mdy(maxt_lapse$date)
    mint_lapse$date <- mdy(mint_lapse$date)        
    maxt_lapse <- maxt_lapse[maxt_lapse$date < ymd("2000-01-01") & maxt_lapse$date > ymd("1969-12-31"),]
    mint_lapse <- mint_lapse[mint_lapse$date < ymd("2000-01-01") & mint_lapse$date > ymd("1969-12-31"),] 
    
    maxt_lapse$day <- yday(maxt_lapse$date)
    mint_lapse$day <- yday(mint_lapse$date)
    
    maxt_lapse$month <- month(maxt_lapse$date)
    mint_lapse$month <- month(mint_lapse$date)
    
    maxt_lapse$year <- year(maxt_lapse$date)
    mint_lapse$year <- year(mint_lapse$date)
    
  # import model coefficients
    maxtmod <- read.csv(paste(climate_path,"micro_mod/maxcoef_out.csv",sep=""),row.names=1)
    mintmod <- read.csv(paste(climate_path,"micro_mod/mincoef_out.csv",sep=""),row.names=1)
    
  # model effect of intercept on lapse rate
     lm_maxt <- lm(slope ~ intercept, maxt_lapse)
     lm_mint <- lm(slope ~ intercept, mint_lapse)
     
     plot(slope ~ intercept,maxt_lapse,ylim=c(-0.022,0.021),xlim=c(-25,42))
     abline(a=lm_maxt$coefficients["(Intercept)"], b=lm_maxt$coefficients["intercept"], col="red")
     
     plot(slope ~ intercept,mint_lapse,ylim=c(-0.022,0.021),xlim=c(-25,42))
     abline(a=lm_mint$coefficients["(Intercept)"], b=lm_mint$coefficients["intercept"], col="red")
     
     maxt_slopeadj <- lm_maxt$coefficients["intercept"] * 4
     mint_slopeadj <- lm_mint$coefficients["intercept"] * 4
     
  # apply model to create +4 degrees scenario
     maxt_lapse_cc <- maxt_lapse
     mint_lapse_cc <- mint_lapse    

     maxt_lapse_cc$intercept <- maxt_lapse_cc$intercept +4
     mint_lapse_cc$intercept <- mint_lapse_cc$intercept +4

     maxt_lapse_cc$slope <- maxt_lapse$slope + maxt_slopeadj
     mint_lapse_cc$slope <- mint_lapse$slope + mint_slopeadj
             
  # import files created by Micro_annualrasters.R           
  
    min_tmp_names <- list.files(paste(climate_path,"tmpMin/",sep=""),full.names=T)
    yearminstack <- stack(min_tmp_names,quick=T)

    
    max_tmp_names <- list.files(paste(climate_path,"tmpMax/",sep=""),full.names=T)
    yearmaxstack <- stack(max_tmp_names,quick=T)

    
  # function to calculate historical and future microclimate from files
    
    CalcDailyClim <- function(year) {
    ## historical 
      ## min temp
      min_lapse <- mint_lapse[mint_lapse$year==year,]
      min_lapse <- min_lapse[order(min_lapse$date),]
      ## max temp
      max_lapse <- maxt_lapse[maxt_lapse$year==year,]
      max_lapse <- max_lapse[order(max_lapse$date),]
      
      
      for(j in 1:365){
        dminlapse <- min_lapse[min_lapse$day == j,]
        
        
        minSYN  <- dminlapse$intercept + (dminlapse$slope * elev)
        
        rad <- radstack[[j]]
        
        min <- yearminstack[[j]] +
          (mintmod["minSYN",1] * minSYN) +
          (mintmod["minSYN:RAD",1] * minSYN * rad) +
          (mintmod["minSYN:ELEV",1] * minSYN * elev) +
          (mintmod["minSYN:LOG.STRDST",1] * minSYN * strdist) +
          (mintmod["minSYN:I(log(TCI))",1] * minSYN * tci)
        
        names(min) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
        
        writeRaster(min,
                    paste(climate_path,"minTs/y_",year,"_d_",sprintf("%03d",j),sep=""),
                    format="GTiff")
        
        
        
        dmaxlapse <- max_lapse[max_lapse$day == j,]
        
        
        maxSYN  <- dmaxlapse$intercept + (dmaxlapse$slope * elev)
        
        
        max <- yearmaxstack[[j]] +
          (maxtmod["maxSYN",1] * maxSYN) +
          (maxtmod["maxSYN:RAD",1] * maxSYN * rad) +
          (maxtmod["maxSYN:ELEV",1] * maxSYN * elev) +
          (maxtmod["maxSYN:TOTRAD",1] * maxSYN * totrad ) +
          (maxtmod["maxSYN:LOG.STRDST",1] * maxSYN * strdist)
        
        names(max) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
        
        writeRaster(max,
                    paste(climate_path,"maxTs/y_",year,"_d_",sprintf("%03d",j),sep=""),
                    format="GTiff")
        
        
        
      }
      
    ## climate change
      ## min temp
      min_lapse <- mint_lapse_cc[mint_lapse_cc$year==year,]
      min_lapse <- min_lapse[order(min_lapse$date),]
      ## max temp
      max_lapse <- maxt_lapse_cc[maxt_lapse_cc$year==year,]
      max_lapse <- max_lapse[order(max_lapse$date),]
      
        
        for(j in 1:365){
          dminlapse <- min_lapse[min_lapse$day == j,]
          
          
          minSYN  <- dminlapse$intercept + (dminlapse$slope * elev)
          
          rad <- radstack[[j]]
          
          min <- yearminstack[[j]] +
            (mintmod["minSYN",1] * minSYN) +
            (mintmod["minSYN:RAD",1] * minSYN * rad) +
            (mintmod["minSYN:ELEV",1] * minSYN * elev) +
            (mintmod["minSYN:LOG.STRDST",1] * minSYN * strdist) +
            (mintmod["minSYN:I(log(TCI))",1] * minSYN * tci)
          
          names(min) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
          
          writeRaster(min,
                      paste(climate_path,"minTs_cc/y_",year,"_d_",sprintf("%03d",j),sep=""),
                      format="GTiff")
          
          
          
          dmaxlapse <- max_lapse[max_lapse$day == j,]
          
          
          maxSYN  <- dmaxlapse$intercept + (dmaxlapse$slope * elev)
          
          
          max <- yearmaxstack[[j]] +
            (maxtmod["maxSYN",1] * maxSYN) +
            (maxtmod["maxSYN:RAD",1] * maxSYN * rad) +
            (maxtmod["maxSYN:ELEV",1] * maxSYN * elev) +
            (maxtmod["maxSYN:TOTRAD",1] * maxSYN * totrad ) +
            (maxtmod["maxSYN:LOG.STRDST",1] * maxSYN * strdist)
          
          names(max) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
          
          writeRaster(max,
                      paste(climate_path,"maxTs_cc/y_",year,"_d_",sprintf("%03d",j),sep=""),
                      format="GTiff")
          
          
          
        }
        
      }
      
    clus <- makePSOCKcluster(15)
    clusterExport(cl=clus,varlist=ls())
    clusterEvalQ(cl=clus,library(raster))
    clusterApply(cl=clus,x=1970:1999,fun=CalcDailyClim)
    stopCluster(clus)
    
  # avoid writng GIS variables back to global env
    allvars <- ls()
    
    rmvars <- allvars[!allvars %in% mainvars]
    
    rm(list=rmvars)
    