####  explore GCM warming compared to 2.5 minute WorldClim ####
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
      
      world_MAT <- raster(paste0(climate_path,"wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif"))

#### crop worldclim data to park, correct for *10 values, and resample to UTM 30m ####
      
  # reduce to park area
    parkbound_tr <- spTransform(park_bound,crs(world_MAT))
    parkMAT <- crop(world_MAT,parkbound_tr)

    
#### check GCMs ####
  GCM_files <- list.files(paste0(climate_path,"GCM/"),full.names=T)
    # ssp explanation here: https://www.carbonbrief.org/cmip6-the-next-generation-of-climate-models-explained
    # data from here: https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html
    # CMIP6 info and terms of use: https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-1.html
    
  GCM_names <- simplify2array(strsplit(GCM_files,"_"))[6,]
    
    
  park_delta <- list()
  
  for(i in 1:length(GCM_files)){
    rast <- stack(GCM_files[i])
    parkrast <- crop(rast[[1]],parkbound_tr)
    park_delta[[i]] <- parkrast - parkMAT
  }
  
  deltaStack <- stack(park_delta)
  names(deltaStack) <- GCM_names
    
  deltaStack_proj <- projectRaster(deltaStack,elev,method="ngb")
  
  deltaStack_proj <- mask(deltaStack_proj,park_bound)
  deltaStack <- mask(deltaStack,parkbound_tr)
  
  plot(deltaStack_proj)
  plot(park_bound,add=T)
  
  col <- colorRampPalette(rev(brewer.pal(9, 'RdGy')))
  
  levelplot(deltaStack_proj,
            maxpixels=1e6,
            margin=F,
            col.regions=col,
            colorkey=list(space="right"),
            at=seq(min(cellStats(deltaStack_proj,"min")),
                   max(cellStats(deltaStack_proj,"max")),
                   len=100),
            scales=list(draw=F)) +
    layer(sp.polygons(park_bound))
  
  
  delta_summary <- data.frame(name=names(deltaStack),
                              mean=cellStats(deltaStack,"mean"),
                              max=cellStats(deltaStack,"max"),
                              min=cellStats(deltaStack,"min"))
  
  delta_summary

  