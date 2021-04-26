#### comparison of model fit across species
## Jordan Stark, Mar 2021


#### setup ####

  # Packages
    library(Hmsc)

  # threshold for inclusion in high_AUC
    thresh <- 0.8

  # paths
    #model_path <- "E:/Smokies_Veg/model_files/"
    #db_path <- "E:/Smokies_Veg/csv_data/"

#### load models of interest and extract AUC by species ####
    load(paste0(model_path,modlist_path))

    calcAUC <- function(mod) {
      pred <- computePredictedValues(mod)
      AUC <- evaluateModelFit(mod,pred)[["AUC"]]
    }
    
    AUC <- data.frame(sapply(modlist,calcAUC))
    AUC <- data.frame(AUC)
    
    AUCdf <- cbind(sp=modlist[[1]]$spNames,AUC)
    names(AUCdf)[2:length(AUCdf)] <- simplify2array(strsplit(names(AUCdf[2:length(AUCdf)]),".",fixed=T))[1,]

    high_AUC <- AUCdf[AUCdf$micro>thresh & AUCdf$fine>thresh & AUCdf$bio>thresh,]
    
  # remove ATHYFIL, predictions are inconsistent and due to positive quadratic term
    high_AUC <- high_AUC[high_AUC$sp != "ATHYFIL",]
    
    write.csv(high_AUC,paste0(output_path,"high_AUC_",short_modname,".csv"))
