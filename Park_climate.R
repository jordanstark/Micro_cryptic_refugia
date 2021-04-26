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
      biorast <- getData(name="worldclim",var="bio",res=0.5,lon=-83,lat=35,path=climate_path)[[1:11]]
      
      park_bound <- readOGR(paste(gis_path,"GRSM_BOUNDARY_POLYGON",sep=""))
      park_bound <- park_bound[park_bound$OBJECTID<18,] # remove external trails etc
      
      elev <- raster(paste0(fridley_GIS,"elev.txt"))
      crs(elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
      
#### crop worldclim data to park, correct for *10 values, and resample to UTM 30m ####
    
  # reduce to park area
    parkbound_tr <- spTransform(park_bound,crs(biorast))
    biopark <- crop(biorast,parkbound_tr)
    #biopark <- mask(biopark,parkbound_tr) this increases the number of NA values around the edges
  
  # isolate values, correct units  
    bio_MAT <- biopark[[1]]/10
    bio_MAXmo <- biopark[[5]]/10
    bio_MINmo <- biopark[[6]]/10
    
      
  # resample with method="ngb" to copy values to 30m grid without splining
    bio_MAT <- projectRaster(bio_MAT,elev,method="ngb")
    bio_MAXmo <- projectRaster(bio_MAXmo,elev,method="ngb")
    bio_MINmo <- projectRaster(bio_MINmo,elev,method="ngb")
    
    
  # save files
    writeRaster(bio_MAT,paste0(climate_path,"bio_MAT"),format="GTiff")
    writeRaster(bio_MAXmo,paste0(climate_path,"bio_MAXmo"),format="GTiff")  
    writeRaster(bio_MINmo,paste0(climate_path,"bio_MINmo"),format="GTiff")
    
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
    
    finebio_MAT <- lapsecalc(bio_MAT)
    finebio_MAXmo <- lapsecalc(bio_MAXmo)
    finebio_MINmo <- lapsecalc(bio_MINmo)
    
  # save files
    writeRaster(finebio_MAT,paste0(climate_path,"fine_MAT"),format="GTiff")
    writeRaster(finebio_MAXmo,paste0(climate_path,"fine_MAXmo"),format="GTiff")
    writeRaster(finebio_MINmo,paste0(climate_path,"fine_MINmo"),format="GTiff")
    
#### calculate future climate bioclim variables
    
  # save MAT, MAXmo, MINmo for +4 degrees
    bio_MAT.cc <- bio_MAT + 4
    bio_MAXmo.cc <- bio_MAXmo + 4
    bio_MINmo.cc <- bio_MINmo + 4
    
    writeRaster(bio_MAT.cc,paste0(climate_path,"bio_MAT_plus4"),format="GTiff")
    writeRaster(bio_MAXmo.cc,paste0(climate_path,"bio_MAXmo_plus4"),format="GTiff")
    writeRaster(bio_MINmo.cc,paste0(climate_path,"bio_MINmo_plus4"),format="GTiff")
    
    fine_MAT.cc <- finebio_MAT + 4
    fine_MAXmo.cc <- finebio_MAXmo + 4
    fine_MINmo.cc <- finebio_MINmo + 4
    
    writeRaster(fine_MAT.cc,paste0(climate_path,"fine_MAT_plus4"),format="GTiff")
    writeRaster(fine_MAXmo.cc,paste0(climate_path,"fine_MAXmo_plus4"),format="GTiff")
    writeRaster(fine_MINmo.cc,paste0(climate_path,"fine_MINmo_plus4"),format="GTiff")

    
#### write parkwide rasters of northing and easting ####
    Easting <- Northing <- elev
      
    values(Easting) <-  rep(xFromCol(Easting),nrow(Easting))   
    values(Northing) <-  rep(yFromRow(Northing),each=ncol(Northing))   
    
    writeRaster(Easting,paste0(gis_path,"Easting"),format="GTiff")
    writeRaster(Northing,paste0(gis_path,"Northing"),format="GTiff")
    