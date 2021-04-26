## prep raster stacks for model input
  library(raster)
  library(rgdal)

# all GIS layers that could be used
# physical GIS data
  elev    <- raster(paste0(fridley_GIS,"elev.txt"))
  tci     <- raster(paste0(fridley_GIS,"tci_cor.asc"))
  totrad  <- raster(paste0(fridley_GIS,"totrad.txt"))
  
  
  crs(elev) <- crs(totrad) <- crs(tci) <- CRS("+proj=utm +zone=17 +datum=NAD27")
  
  log_tci <- log(tci)
  
  names(log_tci) <- "log_tci"
  names(totrad) <- "totrad"
  names(elev) <- "elev"
  
  Northing <- raster(paste0(gis_path,"Northing.tif"))
  Easting <- raster(paste0(gis_path,"Easting.tif"))
  log_PlotSize <- setValues(elev,log(900))
  names(log_PlotSize) <- c("log_PlotSize")
  
  phys_stack <- stack(elev,log_tci,totrad,Northing,Easting,log_PlotSize)

# climate data
  bio_MAT <- raster(paste0(climate_path,"_bio_MAT.tif"))
  bio_MAXmo <- raster(paste0(climate_path,"_bio_MAXmo.tif"))
  bio_MINmo <- raster(paste0(climate_path,"_bio_MINmo.tif"))
  
  fine_MAT <- raster(paste0(climate_path,"_fine_MAT.tif"))
  fine_MAXmo <- raster(paste0(climate_path,"_fine_MAXmo.tif"))
  fine_MINmo <- raster(paste0(climate_path,"_fine_MINmo.tif"))
  
  broad_stack <- stack(bio_MAT,bio_MAXmo,bio_MINmo,
                       fine_MAT,fine_MAXmo,fine_MINmo)
  names(broad_stack) <- c("bio_MAT","bio_MAXmo","bio_MINmo",
                          "fine_MAT","fine_MAXmo","fine_MINmo")
  
  micro_MAT <- raster(paste0(climate_path,"micro_MAT.tif"))
  micro_MAXmo <- raster(paste0(climate_path,"micro_MAXmo.tif"))
  micro_MINmo <- raster(paste0(climate_path,"micro_MINmo.tif"))
  
  # plus4 climate data
  bio_MAT_cc <- raster(paste0(climate_path,"bio_MAT_plus4.tif"))
  bio_MAXmo_cc <- raster(paste0(climate_path,"bio_MAXmo_plus4.tif"))
  bio_MINmo_cc <- raster(paste0(climate_path,"bio_MINmo_plus4.tif"))
  
  fine_MAT_cc <- raster(paste0(climate_path,"fine_MAT_plus4.tif"))
  fine_MAXmo_cc <- raster(paste0(climate_path,"fine_MAXmo_plus4.tif"))
  fine_MINmo_cc <- raster(paste0(climate_path,"fine_MINmo_plus4.tif"))
  
  broad_stack_cc <- stack(bio_MAT_cc,bio_MAXmo_cc,bio_MINmo_cc,
                          fine_MAT_cc,fine_MAXmo_cc,fine_MINmo_cc)
  
  micro_MAT_cc <- raster(paste0(climate_path,"micro_MAT_plus4.tif"))
  micro_MAXmo_cc <- raster(paste0(climate_path,"micro_MAXmo_plus4.tif"))
  micro_MINmo_cc <- raster(paste0(climate_path,"micro_MINmo_plus4.tif"))
  

# crop bio and fine to same ext as micro data
  ext <- extent(micro_MAT)
  
  broad_stack <- crop(broad_stack,ext)
  broad_stack_cc <- crop(broad_stack_cc,ext)
  
  phys_stack <- crop(phys_stack,ext)
  
  
  all_stack <- stack(phys_stack,
                     broad_stack,
                     micro_MAT,micro_MAXmo,micro_MINmo)
  all_stack_cc <- stack(phys_stack,
                        broad_stack_cc,
                        micro_MAT_cc,micro_MAXmo_cc,micro_MINmo_cc)
  names(all_stack_cc) <- gsub("plus4",names(all_stack_cc),replacement="cc")
  
# crop all to park boundary
  parkbound <- readOGR(paste(gis_path,"GRSM_BOUNDARY_POLYGON",sep=""))
  parkbound <- spTransform(parkbound,crs(all_stack))
  parkbound <- parkbound[parkbound$OBJECTID<18,]

  all_stack <- mask(all_stack,parkbound)
  all_stack_cc <- mask(all_stack_cc,parkbound)
      
# save data
  if(!dir.exists(paste0(model_inputs,"model_GIS/"))) dir.create(paste0(model_inputs,"model_GIS/"))
  if(!dir.exists(paste0(model_inputs,"model_GIS_cc/"))) dir.create(paste0(model_inputs,"model_GIS_cc/"))
  
  
  writeRaster(all_stack,
              filename=paste0(model_inputs,"model_GIS/"),
              format="GTiff",
              bylayer=T, suffix="names")
  writeRaster(all_stack_cc,
              filename=paste0(model_inputs,"model_GIS_cc/"),
              format="GTiff",
              bylayer=T, suffix="names")
  