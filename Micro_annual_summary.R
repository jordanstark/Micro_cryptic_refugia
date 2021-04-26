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
  
                  
#### calculate MAT, MAXmo, MINmo for each year and average ####
    maxt_files <- list.files(paste(climate_path,"maxTs/",sep=""),full.names=T)
    mint_files <- list.files(paste(climate_path,"minTs/",sep=""),full.names=T)
    
    
    maxt_cc_files <- list.files(paste(climate_path,"maxTs_cc/",sep=""),full.names=T)
    mint_cc_files <- list.files(paste(climate_path,"minTs_cc/",sep=""),full.names=T)
    
    MonthAssg <- month(as_date(0:364))
    
    
    calcBioVars <- function(year) {
      
        year_max <- grep(year,maxt_files)
        year_min <- grep(year,mint_files)
        
        maxt <- stack(maxt_files[year_max],quick=T)
        mint <- stack(mint_files[year_min],quick=T)
        
        avgmax <- vector(mode="list",length=12)
        avgmin <- vector(mode="list",length=12)
        
        for(j in 1:12) {
          max_mo <- maxt[[which(MonthAssg == j)]]
          min_mo <- mint[[which(MonthAssg == j)]]
          
          avgmax[[j]] <- calc(max_mo,mean)
          avgmin[[j]] <- calc(min_mo,mean)
        }
        
        avgmax <- stack(avgmax)
        avgmin <- stack(avgmin)
        
        
        # calculate bio1 (mean temp), bio5 (max temp warmest mo) and bio6 (min temp warmest mo)
        tavg <- (avgmax + avgmin) / 2
        
        MAT <- calc(tavg, mean, filename=paste0(climate_path,"summary/MAT_y_",year,".tif"))
        MAXmo <- calc(avgmax,max,filename=paste0(climate_path,"summary/MAXmo_y_",year,".tif"))
        MINmo <- calc(avgmin,min,filename=paste0(climate_path,"summary/MINmo_y_",year,".tif"))
        
        
        ## calculate for +4 degrees scenario
        year_max_cc <- grep(year,maxt_cc_files)
        year_min_cc <- grep(year,mint_cc_files)
        
        maxt_cc <- stack(maxt_cc_files[year_max_cc],quick=T)
        mint_cc <- stack(mint_cc_files[year_min_cc],quick=T)
        
        avgmax_cc <- vector(mode="list",length=12)
        avgmin_cc <- vector(mode="list",length=12)
        
        for(j in 1:12) {
          max_mo <- maxt_cc[[which(MonthAssg == j)]]
          min_mo <- mint_cc[[which(MonthAssg == j)]]
          
          avgmax_cc[[j]] <- calc(max_mo,mean)
          avgmin_cc[[j]] <- calc(min_mo,mean)
        }
        
        avgmax_cc <- stack(avgmax_cc)
        avgmin_cc <- stack(avgmin_cc)
        
        
        # calculate MAT and SEAS following biovars() #bio1, bio2 and bio4
        tavg_cc <- (avgmax_cc + avgmin_cc) / 2
        
        MAT_cc <- calc(tavg_cc, mean, filename=paste0(climate_path,"summary_cc/MAT_y_",year,".tif"))
        MAXmo_cc <- calc(avgmax_cc,max,filename=paste0(climate_path,"summary_cc/MAXmo_y_",year,".tif"))
        MINmo_cc <- calc(avgmin_cc,min,filename=paste0(climate_path,"summary_cc/MINmo_y_",year,".tif"))
        
    }
    
    
    clus <- makePSOCKcluster(15)
    clusterExport(cl=clus,varlist=ls())
    clusterEvalQ(cl=clus,library(raster))
    clusterApply(cl=clus,x=1970:1999,fun=calcBioVars)
    stopCluster(clus)
    
    # avoid writng GIS variables back to global env
    allvars <- ls()
    
    rmvars <- allvars[!allvars %in% mainvars]
    
    rm(list=rmvars)
    
    
    