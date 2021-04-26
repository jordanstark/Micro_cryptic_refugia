#### Prep database files for species distribution analysis
## Jordan Stark
## Winter 2020-2021

#### Setup ####
  # packages
    library(tidyr) # for pivot_longer and _wider
    library(raster) # raster data
    library(sp) # spatial points
    library(rgdal)




#### import data ####
  # database files
    sp_records      <- read.csv(paste(db_path,"HerbdataSVD_3-14.csv",sep=""))
    sp_records$ID   <- NULL
    sp_info         <- read.csv(paste(db_path,"CarSpList_SVD.csv",sep=""))
    sp_info$ID      <- NULL
    plot_records    <- read.csv(paste(db_path,"Rasters_May21.csv",sep=""))
    plot_records$ID <- NULL

      
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
    
    phys_stack <- stack(elev,log_tci,totrad,Northing,Easting)
  
    
  # climate data
    bio_MAT <- raster(paste0(climate_path,"_bio_MAT.tif"))
    
    fine_MAT <- raster(paste0(climate_path,"_fine_MAT.tif"))
    
    micro_MAT <- raster(paste0(climate_path,"micro_MAT.tif"))

    ext <- extent(micro_MAT)
    
  # crop bio and fine to same ext as micro data
    bio_MAT <- crop(bio_MAT,ext)
    
    fine_MAT <- crop(fine_MAT,ext)
    
    clim_stack <- stack(bio_MAT,fine_MAT,micro_MAT)
    names(clim_stack) <- c("bio_MAT","fine_MAT","micro_MAT")

#### combine & clean database ####

  # create combined df with all sp info and plot data
    info_good <- sp_info[,c("SppID","FAMILY","NC_CODE","US_COMMON_NAME","NC_Std")]
    info_good$NC_Std <- as.numeric(info_good$NC_Std)
    
    all_db <- merge(sp_records,info_good,all.x=T)

    
  # remove species with ambiguous IDs
    all_db <- all_db[all_db$NC_Std %in% c(1,5,6),]
    
    all_db$NC_Std <- NULL
    
  # check number of records per sp
    numrec <- aggregate(coverCVS ~ NC_CODE, all_db,length)
    names(numrec) <- c("NC_CODE","recs")
    
    raresp <- numrec$NC_CODE[numrec$recs < 50]
    
  # remove rare species (<10 occurences)
    all_db <- all_db[-which(all_db$NC_CODE %in% raresp),]
    
    
    
#### transform to occurence records ####
    
  # spread cover categories to wide format
    sp_IDs <- c("SppID","GENUS","GEN.SPP","FAMILY","US_COMMON_NAME")
      # need to remove this info to avoid duplicates
    only_ID <- all_db[,!names(all_db) %in% sp_IDs]
    
    
    
    bad_rows <- which(duplicated(only_ID[,c("PLOT_ID","Date","NC_CODE")]))
    only_ID <- only_ID[-bad_rows,] #removes AGERALTA listed twice in plot 
    
    wide <- pivot_wider(only_ID,names_from=NC_CODE,
                        values_from=coverCVS,values_fill=0)
    
  # set zeros to 'NA' in surveys that were type 'T' (trees only)
    
    for(i in 1:length(wide[1,])){
      ST <- wide$Stratum.Type[i]
      
      if(ST=="T"|ST=="T ") {
        wide[i,][wide[i,]==0] <- NA
      }
    }
    
  # remove earlier surveys of plots that were re-measured later
    wide$year <- as.numeric(simplify2array(strsplit(wide$Date,split="-",fixed=T))[3,])
    wide <- wide[order(wide$year),]
    repeated_plots <- which(duplicated(wide$PLOT_ID,fromLast=T))
    
    wide <- wide[-repeated_plots,]
    
    veg_plots <- data.frame(wide[,c("PLOT_ID","Date","PlotSize")])
    
    wide <- wide[,!names(wide) %in% c("year","PlotSize","Stratum.Type")]
    
    
#### extract all variables which will be used in Hmsc models ####
    
  # make plot records into spatial points
    veg_plots$log_PlotSize <- log(veg_plots$PlotSize)
    veg_plots$PlotSize <- NULL
    names(veg_plots) <- c("PLOT_ID","DATE","log_PlotSize")
    
    veg_plot_recs <- merge(veg_plots,
                           plot_records[,c("PLOT_ID","DATE","UTM_E","UTM_N")],
                           sort=F,
                           all.x=T)
    spatial_plots <- SpatialPointsDataFrame(coords = veg_plot_recs[,c("UTM_E","UTM_N")],
                                            proj4string=CRS("+proj=utm +zone=17 +datum=NAD27"),
                                            data=veg_plot_recs[,c("PLOT_ID","UTM_E","UTM_N","log_PlotSize")])

    
  # Extract GIS values at plots
    phys_points <- raster::extract(phys_stack,spatial_plots)
    clim_points <- raster::extract(clim_stack,spatial_plots)
    
    point_dat <- data.frame(cbind(phys_points,clim_points))
    point_dat <- cbind(spatial_plots@data,point_dat)

    point_dat <- unique(point_dat)
    
  # remove points with NA values
    
    complete_plots <- complete.cases(point_dat)
    
    point_dat <- point_dat[complete_plots,]
    abund_mat <- wide[complete_plots,] #remove bad plots from veg data too
    
  # create model matrix with all possible terms

    xvars <- point_dat[5:length(point_dat)]
    xvars.sq <- xvars^2
    names(xvars.sq) <- paste0(names(xvars),"2")
    
    allx_noint <- cbind(xvars,xvars.sq)
    
    allx <- data.frame(model.matrix(~.+.^2,data=allx_noint))
    allx$log_PlotSize <- point_dat$log_PlotSize
    
    allx_plots <- cbind(point_dat[,1:3],allx)
    
  # make sure this is in the same order as the veg data
    stopifnot(identical(allx_plots$PLOT_ID,abund_mat$PLOT_ID))
    
  # center and scale all variables
    scalemeans <- apply(allx,2,mean)
    scalesd    <- apply(allx,2,sd)
    
    scaled_mat <- data.frame(matrix(nrow=length(allx[,1]),ncol=length(allx)))
    names(scaled_mat) <- names(allx)
    
    for(i in 1:length(scaled_mat[,1])) {
      scaled_mat[i,] <- (allx[i,] - scalemeans) / scalesd
    }
    
    # save scale means and sd
    scales <- data.frame(means=scalemeans,
                         sds=scalesd)
    scales$par <- row.names(scales)
    

  
  # transform veg data to presence/absence
    presence_mat <- abund_mat[3:length(abund_mat)]
    presence_mat[presence_mat>0] <- 1
    presence_mat <- cbind(abund_mat[,1:2],presence_mat)
 

#### save data ####
  write.csv(presence_mat,paste0(model_inputs,"presence_data.csv"),row.names=F) 
  write.csv(scaled_mat,paste0(model_inputs,"env_scaled.csv"),row.names=F)
  write.csv(scales,paste0(model_inputs,"scalepars.csv"),row.names=F)
   