#### prep spatial data across GSMNP for veg analysis ####
#### setup ####
  # packages
    library(raster)
    library(sp)
    library(rgdal)



  # data paths
      #climate_path <- "E:/Smokies_Veg/climate_data/"  
      #gis_path     <- "E:/Smokies_Veg/phys_GIS/" 
      #fridley_GIS <- "E:/Smokies_Veg/phys_GIS/gsmnp_ascii/"
      
  # import climate and park data
      park_bound <- readOGR(paste(gis_path,"GRSM_BOUNDARY_POLYGON",sep=""))
      park_bound <- park_bound[park_bound$OBJECTID<18,] # remove external trails etc
      
      elev <- raster(paste0(fridley_GIS,"elev.txt"))
      crs(elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
      
      world_MAT <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
      world_DIU <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_2.tif"))
      world_ISO <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_3.tif"))
      world_SEAS <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_4.tif"))
      world_MAXmo <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_5.tif"))
      world_MINmo <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_6.tif"))
      world_ANR <- raster(paste0(climate_path,"wc2.1_30s_bio/wc2.1_30s_bio_7.tif"))
      
      
      world_stack <- stack(world_MAT,world_DIU,world_ISO,world_SEAS,world_MAXmo,world_MINmo,world_ANR)
      names(world_stack) <- c("bio_MAT","bio_DIU","bio_ISO","bio_SEAS","bio_MAXmo","bio_MINmo","bio_ANR")      
      
#### crop worldclim data to park, correct for *10 values, and resample to UTM 30m ####
      
  # reduce to park area
    parkbound_tr <- spTransform(park_bound,crs(world_stack))
    park_stack <- crop(world_stack,parkbound_tr)

      
  # resample with method="ngb" to copy values to 30m grid without splining
    park_stack <- projectRaster(park_stack,elev,method="ngb")
    
    
  # save files
    writeRaster(park_stack,bylayer=T,filename=paste0(climate_path),format="GTiff",suffix="names")
    
    
#### calculate lapse-downscaled bioclim variables ####
    
  bioclim_elev <- raster(paste(climate_path,"wc2.1_30s_elev/wc2.1_30s_elev.tif",sep=""))
  bioclim_elev_park <- projectRaster(bioclim_elev,elev,method="ngb")
  
  elev.diff <- bioclim_elev_park - elev
    
  # calculate correlation with elevation and apply to finescale elevation

    lapsecalc <- function(clim) {
      rasts <- stack(clim,bioclim_elev_park)
      df <- data.frame(values(rasts))
      names(df) <- c("clim","elev")
      
      mod <-  lm(clim ~ elev, data=df)
      
      slope <- mod$coefficients["elev"]
      int <- mod$coefficients["(Intercept)"]
    
      clim_ds <- clim - (slope*elev.diff) 
      
    }
    
    fine_stack <- list()
    
    for(i in 1:nlayers(park_stack)){
      fine_stack[[i]] <- lapsecalc(park_stack[[i]])
    }
    
    fine_stack <- stack(fine_stack)
    names(fine_stack) <- sub("bio","fine",names(park_stack))
    
  # save files
    writeRaster(fine_stack,bylayer=T,filename=paste0(climate_path),format="GTiff",suffix="names")
    

#### calculate future climate bioclim variables
    
  # save MAT, MAXmo, MINmo for +4 degrees
    bio_MAT.cc <- park_stack$bio_MAT + 4
    bio_MAXmo.cc <- park_stack$bio_MAXmo + 4
    bio_MINmo.cc <- park_stack$bio_MINmo + 4
    
    writeRaster(bio_MAT.cc,paste0(climate_path,"bio_MAT_plus4"),format="GTiff")
    writeRaster(bio_MAXmo.cc,paste0(climate_path,"bio_MAXmo_plus4"),format="GTiff")
    writeRaster(bio_MINmo.cc,paste0(climate_path,"bio_MINmo_plus4"),format="GTiff")
    
    fine_MAT.cc <- fine_stack$fine_MAT + 4
    fine_MAXmo.cc <- fine_stack$fine_MAXmo + 4
    fine_MINmo.cc <- fine_stack$fine_MINmo + 4
    
    writeRaster(fine_MAT.cc,paste0(climate_path,"fine_MAT_plus4"),format="GTiff")
    writeRaster(fine_MAXmo.cc,paste0(climate_path,"fine_MAXmo_plus4"),format="GTiff")
    writeRaster(fine_MINmo.cc,paste0(climate_path,"fine_MINmo_plus4"),format="GTiff")

    
#### write parkwide rasters of northing and easting ####
    Easting <- Northing <- elev
      
    values(Easting) <-  rep(xFromCol(Easting),nrow(Easting))   
    values(Northing) <-  rep(yFromRow(Northing),each=ncol(Northing))   
    
    writeRaster(Easting,paste0(gis_path,"Easting"),format="GTiff")
    writeRaster(Northing,paste0(gis_path,"Northing"),format="GTiff")
    