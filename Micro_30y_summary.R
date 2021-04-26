#### run Fridley microclimate model; adjust climate predictions for change in lapse

#### setup ####
  # packages
    library(raster)
    library(rgdal)

  # paths
    #fridley_GIS <- "E:/Smokies_Veg/phys_GIS/gsmnp_ascii/"
    #climate_path <- "E:/Smokies_Veg/climate_data/" 
  
      
#### calculate mean MAT, MAXmo, MINmo, over 30 y ####

  calc30mean <- function(var,filelist,out_name){
    annual <- stack(grep(var,filelist,value=T),quick=T)
    
    avg30 <- calc(annual,mean)
    
    filled <- focal(avg30, w=matrix(1,3,3), fun=mean, NAonly=T,na.rm=T)
    names(filled) <- var
    
    writeRaster(filled,paste0(climate_path,out_name,".tif"))
  }

    
  # historical climate
    summary.files <- list.files(paste(climate_path,"summary/",sep=""),full.names=T)
    
    calc30mean("MAT",summary.files,"micro_MAT")
    calc30mean("MAXmo",summary.files,"micro_MAXmo")
    calc30mean("MINmo",summary.files,"micro_MINmo")
    
  # +4 degrees
    cc.summary.files <- list.files(paste(climate_path,"summary_cc/",sep=""),full.names=T)
  
    calc30mean("MAT",cc.summary.files,"micro_MAT_plus4")
    calc30mean("MAXmo",cc.summary.files,"micro_MAXmo_plus4")
    calc30mean("MINmo",cc.summary.files,"micro_MINmo_plus4")
  