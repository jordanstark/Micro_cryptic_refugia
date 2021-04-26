#### apply coefficients of HMSC models to predict distribution across park
## and under plus 4 degrees scenario



#### setup ####
  # packages
    library(Hmsc)
    library(sp)
    library(raster)
    library(rgdal)
    library(parallel)


  # directories
    # intermediate_path <- "E:/Smokies_Veg/model_runs/"
    # climate_path <- "E:/Smokies_Veg/climate_data/" 
    # fridley_GIS <- "E:/Fridley_GIS/gsmnp_ascii/"
    # output_path <- "E:/Smokies_Veg/outputs/"

    predicts_path <- paste0(output_path,"predicts_",short_modname,"/")
    if(!dir.exists(predicts_path)) dir.create(predicts_path)

  # information about the models to compare
    load(paste0(model_path,modlist_path))
    
    
    scalepars <- read.csv(paste0(model_inputs,"scalepars.csv"))  

    scaleenv <- function(raster,param) {
      rev.param <- ifelse(length(strsplit(param,".",fixed=T)[[1]])>1,
                          paste0(strsplit(param,".",fixed=T)[[1]][2],
                                 ".",
                                 strsplit(param,".",fixed=T)[[1]][1]),
                          param)

      pars <- scalepars[scalepars$par==param | scalepars$par==rev.param, ]
      
      scaled_raster <- (raster - pars$mean) / pars$sd
    }
    
    spnames <- read.csv(paste0(output_path,"high_AUC_",short_modname,".csv"))[,2]

  # parkwide rasters of model data
    now.files <- list.files(paste0(model_inputs,"model_GIS/"),
                            full.names=T)
    cc.files <- list.files(paste0(model_inputs,"model_GIS_cc/"),
                           full.names=T)
    all_stack <- stack(now.files)
    names(all_stack) <- gsub("X_",names(all_stack),replacement="")
    all_stack_cc <- stack(cc.files)
    names(all_stack_cc) <- gsub("X_",names(all_stack_cc),replacement="")
    names(all_stack_cc) <- gsub("_cc",names(all_stack_cc),replacement="")
    

#### calculate values needed in each model and make predictions

    
  MakeStack <- function(model,scenario) {
    model.params <- attr(terms(model$XFormula),"term.labels")
    
    model_stack <- vector(mode="list",length=length(model.params))
    
    if(scenario=="now") {
      full_stack <- all_stack
    } else if(scenario=="cc") {
      full_stack <- all_stack_cc
    } else stop("check scenario names")
    
    
    for(j in 1:length(model.params)){
      if(model.params[j] %in% names(full_stack)){
        model_stack[[j]] <- full_stack[[model.params[j]]]
      } else if(model.params[j] %in% paste0(names(full_stack),"2")){
        non.sq.name <- strsplit(model.params[j],"2")[[1]][1]
        model_stack[[j]] <- full_stack[[non.sq.name]]^2
      } else {
        layer1_name <- strsplit(model.params[j],".",fixed=T)[[1]][1]
        layer2_name <- strsplit(model.params[j],".",fixed=T)[[1]][2]
        if(layer1_name %in% names(full_stack)){
          layer1 <- full_stack[[layer1_name]]
        } else if(layer1_name %in% paste0(names(full_stack),"2")){
          non.sq.name <- strsplit(layer1_name,"2")[[1]][1]
          layer1 <- full_stack[[non.sq.name]]^2
        } else  stop("error:layers not found, check names")
        if(layer2_name %in% names(full_stack)){
          layer2 <- full_stack[[layer2_name]]
        } else if(layer2_name %in% paste0(names(full_stack),"2")){
          non.sq.name <- strsplit(layer2_name,"2")[[1]][1]
          layer2 <- full_stack[[non.sq.name]]^2
        } else  stop("error:layers not found, check names")
        
        model_stack[[j]] <- layer1 * layer2
      } 
      
      model_stack[[j]] <- scaleenv(model_stack[[j]],model.params[j])
    }
    
    model_stack <- stack(model_stack)
    names(model_stack) <- model.params
    model_stack
  }
  
  model_stack_list <- lapply(modlist,MakeStack,scenario="now")
  model_stack_list_cc <- lapply(modlist,MakeStack,scenario="cc")
  
  all_stack <- NULL
  all_stack_cc <- NULL
  
  extractMed <- function(mod){
    df <- data.frame(t(getPostEstimate(mod,"Beta",q=0.5)$q))
    names(df) <- sub(":",".",mod$covNames,fixed=T)
    df <- cbind(sp=mod$spNames,df)
    df
  }

  parlist <- lapply(modlist,extractMed)
  
  PredictSp <- function(sp,modname,scenario){
    model <- modlist[[modname]]
    pars <- parlist[[modname]]
    if(scenario=="now"){
      modstack <- model_stack_list[[modname]]
    }
    if(scenario=="cc"){
      modstack <- model_stack_list_cc[[modname]]
    }
  
    
    if(identical(names(modstack),names(pars)[3:length(pars)]) !=T){
      stop("parameters and raster layers do not match")
    }
              

    spval <- as.numeric(pars[pars$sp==sp,3:length(pars)])
      
      
    out <- pars[pars$sp==sp,"(Intercept)"]

    for(j in 1:length(spval)) {
      out <- out + spval[j] * modstack[[j]]
    }

    # probit models are in Z-scale, need to transform to probability
    if(model$distr[1,"family"]==2){
      out <- calc(out,pnorm)
    }
      
    writeRaster(out,
                paste0(predicts_path,
                       modname,".",sp,".",scenario,".tif"))
    
    
  }
      
  ## calculate predictions, in parallel
  gc()
  

  clus <- makePSOCKcluster(3)
    #using fewer nodes due to high mem reqs

  clusterEvalQ(cl=clus,library(raster))
  clusterExport(cl=clus,varlist=c("modlist","parlist",
                                  "model_stack_list","model_stack_list_cc",
                                  "predicts_path"))

  for(i in 1:length(modlist)){
    clusterExport(cl=clus,varlist=c("i"))
    clusterApply(cl=clus,x=spnames,
                 fun=PredictSp,modname=names(modlist)[i],
                 scenario="now")
    clusterApply(cl=clus,x=spnames,
                 fun=PredictSp,modname=names(modlist)[i],
                 scenario="cc")
  }
  
  stopCluster(clus)
  
  # correct distribution of ABIEFRA, which is incorrectly predicted a low elevations under climate change 
  # because it has no upper-elevation limit in the park
  
  ABIEFRA_stack <- stack(list.files(predicts_path,pattern="ABIEFRA",full.names=T))
  bio.cc <- ABIEFRA_stack[[grep("bio.*.cc",names(ABIEFRA_stack))]]
  bio.now <- ABIEFRA_stack[[grep("bio.*.now",names(ABIEFRA_stack))]]
  fine.cc <- ABIEFRA_stack[[grep("fine.*.cc",names(ABIEFRA_stack))]]
  fine.now <- ABIEFRA_stack[[grep("fine.*.now",names(ABIEFRA_stack))]]
  micro.cc <- ABIEFRA_stack[[grep("micro.*.cc",names(ABIEFRA_stack))]]
  micro.now <- ABIEFRA_stack[[grep("micro.*.now",names(ABIEFRA_stack))]]
  
  bio.cc.cor <- mask(bio.cc,bio.cc>bio.now,maskvalue=1,updatevalue=0)
  fine.cc.cor <- mask(fine.cc,fine.cc>fine.now,maskvalue=1,updatevalue=0)
  micro.cc.cor <- mask(micro.cc,micro.cc>micro.now,maskvalue=1,updatevalue=0)
  
  writeRaster(bio.cc.cor,paste0(predicts_path,"bio.pa.ABIEFRA.cc.tif"),overwrite=T)
  writeRaster(fine.cc.cor,paste0(predicts_path,"fine.pa.ABIEFRA.cc.tif"),overwrite=T)
  writeRaster(micro.cc.cor,paste0(predicts_path,"micro.pa.ABIEFRA.cc.tif"),overwrite=T)
  
  
  # avoid writng GIS variables back to global env
  allvars <- ls()
  
  rmvars <- allvars[!allvars %in% mainvars]
  
  rm(list=rmvars)
  