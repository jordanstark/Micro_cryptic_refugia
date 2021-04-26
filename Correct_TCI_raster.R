#### change 'NA' TCI areas in tci.txt to 95th percentile


#### setup ####
  # packages
    library(raster)

  # paths
    #fridley_GIS


#### import TCI and remove '25' values (NAs) ####

  tci <- raster(paste0(fridley_GIS,"tci.txt"))
  plot(tci)
  

  # set NA values
  tci_raw <- tci
  tci_raw[tci_raw==25] <- NA  
  plot(tci_raw)  
  
  # calculate 95th percentile in park
  tci.95 <- quantile(tci_raw,probs=0.95)
  
  # apply 95 percentile to fill NAs
  tci_cor <- tci
  tci_cor[tci_cor==25] <- tci.95  
  plot(tci_cor)  
  
  
#### save raster ####
  writeRaster(tci_cor,paste0(fridley_GIS,"tci_cor.txt"),format="ascii")
  