## summarize current and climate change distributions

#### setup ####
  # packages
    library(raster)
    library(rgdal)

  # paths
    #output_path
    predicts_path <- paste0(output_path,"predicts_",short_modname,"/")
    summary_path <- paste0(output_path,"spsummary_",short_modname,"/")
    stable_path <- paste0(output_path,"stable_",short_modname,"/")
    
    if(!dir.exists(summary_path)) dir.create(summary_path)
    if(!dir.exists(stable_path)) dir.create(stable_path)
    
  # setup data files
    predicts <- list.files(paste0(predicts_path),
                           full.names=T)

    now.files <- grep("now",predicts,value=T)
    cc.files <- grep("cc",predicts,value=T)


  # extract sp names
    sp <- unique(simplify2array(strsplit(now.files,split=".",fixed=T))[3,])

#### calculate total expected occupied area for each sp ####
  calcArea <- function(raster){
    cellStats(raster,stat="sum",na.rm=T) * 900 / 1000000
  }

  now.area <- data.frame(sp=sp,
                         micro=NA,
                         fine=NA,
                         bio=NA)

  for(i in 1:length(now.area[,1])){
    sp_i <- now.area$sp[i]
    now.area$micro[i] <- calcArea(raster(paste0(predicts_path,
                                                "micro.pa.",
                                                sp_i,
                                                ".now.tif")))
    now.area$fine[i] <- calcArea(raster(paste0(predicts_path,
                                                "fine.pa.",
                                                sp_i,
                                                ".now.tif")))
    now.area$bio[i] <- calcArea(raster(paste0(predicts_path,
                                                "bio.pa.",
                                                sp_i,
                                                ".now.tif")))
  }

  cc.area <- data.frame(sp=sp,
                         micro=NA,
                         fine=NA,
                         bio=NA)

  for(i in 1:length(cc.area[,1])){
    sp_i <- now.area$sp[i]
    cc.area$micro[i] <- calcArea(raster(paste0(predicts_path,
                                                "micro.pa.",
                                                sp_i,
                                                ".cc.tif")))
    cc.area$fine[i] <- calcArea(raster(paste0(predicts_path,
                                               "fine.pa.",
                                               sp_i,
                                               ".cc.tif")))
    cc.area$bio[i] <- calcArea(raster(paste0(predicts_path,
                                              "bio.pa.",
                                              sp_i,
                                              ".cc.tif")))
  }


  change_area <- data.frame(sp=sp,
                            micro=cc.area$micro - now.area$micro,
                            fine=cc.area$fine - now.area$fine,
                            bio=cc.area$bio - now.area$bio)


  write.csv(now.area,paste0(output_path,"now_area_",short_modname,".csv"))
  write.csv(cc.area,paste0(output_path,"cc_area_",short_modname,".csv"))
  write.csv(change_area,paste0(output_path,"change_area_",short_modname,".csv"))



  calcNumSp <- function(scale,scenario){
    if(scenario=="now"){
      filelist <- now.files
    } else if(scenario=="cc"){
      filelist <- cc.files
    }

    calc(stack(grep(scale,filelist,value=T)),
         fun=sum,na.rm=T,
         filename=paste0(summary_path,
                         "num.",scale,".",scenario,".tif"))
  }

  Num.micro.now <- calcNumSp("micro","now")
  Num.fine.now  <- calcNumSp("fine","now")
  Num.bio.now   <- calcNumSp("bio","now")

  Num.micro.cc <- calcNumSp("micro","cc")
  Num.fine.cc  <- calcNumSp("fine","cc")
  Num.bio.cc   <- calcNumSp("bio","cc")

  Num.micro.change <- Num.micro.cc - Num.micro.now
  Num.fine.change <- Num.fine.cc - Num.fine.now
  Num.bio.change <- Num.bio.cc - Num.bio.now

  writeRaster(Num.micro.change,paste0(summary_path,
                                      "num.micro.change",".tif"))
  writeRaster(Num.fine.change,paste0(summary_path,
                                      "num.fine.change",".tif"))
  writeRaster(Num.bio.change,paste0(summary_path,
                                      "num.bio.change",".tif"))

  stable.area <- data.frame(sp=sp,
                            micro=NA,
                            fine=NA,
                            bio=NA)

  for(i in 1:length(stable.area[,1])){
    sp_i <- stable.area$sp[i]
    micro_part_path <- paste0(predicts_path,"micro.pa.",sp_i)
    fine_part_path <- paste0(predicts_path,"fine.pa.",sp_i)
    bio_part_path <- paste0(predicts_path,"bio.pa.",sp_i)

    micro <- min(raster(paste0(micro_part_path,".now.tif")),
                 raster(paste0(micro_part_path,".cc.tif")))
    writeRaster(micro,filename= paste0(stable_path,"micro_",sp_i,".tif"))

    fine <- min(raster(paste0(fine_part_path,".now.tif")),
                raster(paste0(fine_part_path,".cc.tif")))
    writeRaster(fine,filename= paste0(stable_path,"fine_",sp_i,".tif"))

    bio <- min(raster(paste0(bio_part_path,".now.tif")),
               raster(paste0(bio_part_path,".cc.tif")))
    writeRaster(bio,filename= paste0(stable_path,"bio_",sp_i,".tif"))


    stable.area$micro[i] <- calcArea(micro)
    stable.area$fine[i] <- calcArea(fine)
    stable.area$bio[i] <- calcArea(bio)

  }

  stable.area.frac <- data.frame(sp=stable.area$sp,
                                 micro=stable.area$micro / now.area$micro,
                                 fine=stable.area$fine / now.area$fine,
                                 bio=stable.area$bio / now.area$bio)


  write.csv(stable.area,paste0(output_path,"stable_area_",short_modname,".csv"))
  write.csv(stable.area.frac,paste0(output_path,"stable_area_frac",short_modname,".csv"))


  calcStableNumSp <- function(scale){
    filelist <- list.files(stable_path,full.names=T)
    sum_stable <- calc(stack(grep(scale,filelist,value=T)),
                       fun=sum,na.rm=T)
    sum_now <- raster(paste0(summary_path,"num.",scale,".now.tif"))

    frac_stable <- overlay(sum_stable,sum_now,fun=function(a,b)a/b,
                           filename=paste0(summary_path,"frac.",scale,".stable.tif"))

  }

  micro.stable <- calcStableNumSp("micro")
  fine.stable <- calcStableNumSp("fine")
  bio.stable <- calcStableNumSp("bio")


