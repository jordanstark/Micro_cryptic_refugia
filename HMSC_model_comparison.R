#### evaluation of model fit and predictions from 'Clean_HMSC_models.R'



#### setup ####
  # packages
    library(Hmsc)
    library(ggplot2)
    library(patchwork)

  # directories
    #model_path <- "E:/Smokies_Veg/model_runs/"
    #climate_path <- "E:/Smokies_Veg/climate_data/" 
    #fridley_GIS <- "E:/Smokies_Veg/phys_GIS/gsmnp_ascii/"
    #output_path <- "E:/Smokies_Veg/outputs/"



#### function to analyze and compare models ####
  
  CompareMods <- function(filename){
    
    load(paste0(model_path,filename))

  # check convergence
    
    CheckConv <- function(model){
      m.coda <- convertToCodaObject(model,Gamma=F,Psi=F,Delta=F)
      
      psrf.beta <- fivenum(gelman.diag(m.coda$Beta,multivariate=F)$psrf)

      names(psrf.beta) <- c("Min","Q1","Med","Q3","Max")
      
      return(psrf.beta)
      
    }
    
    conv <- sapply(modlist,CheckConv)    
    conv
    
    max_psrf <- conv["Max",]

    
  # compare WAIC

    WAIC <- sapply(modlist,computeWAIC)
    WAIC
    
  # compare AUC for pa mods
    calcAUC <- function(mod) {
      pred <- computePredictedValues(mod)
      AUC <- evaluateModelFit(mod,pred)[["AUC"]]
    }
    
    AUC <- sapply(modlist,calcAUC)
    
    meanAUC <- data.frame(AUC=colMeans(AUC))
    
    meanAUC$mod <- row.names(meanAUC)


      
  # model comparison
    mod_comp <- data.frame(max_psrf,WAIC)
    mod_comp$mod <- row.names(mod_comp)
    mod_comp <- merge(mod_comp,meanAUC)

    mod_comp$scale <- simplify2array(strsplit(mod_comp$mod,"\\."))[1,]
    mod_comp$formula <- as.character(modlist[[1]]$XFormula[2])
    
    mod_comp
  }

    
#### apply function to each model ####
  
  mod_comp_list <- vector(mode="list",length=length(models))
  
    
  for(i in 1:length(mod_comp_list)) {
    mod_comp_list[[i]] <- CompareMods(models[[i]])
    mod_comp_list[[i]]$desc <- names(models)[i]
  }
    
  mod_comps <- do.call(rbind,mod_comp_list)
  
  mod_comps$desc <- factor(mod_comps$desc,
                           levels=names(models),
                           ordered=T)

  write.csv(mod_comps,paste0(output_path,"model_comparisons.csv"))

#### compare results across model versions ####
  WAIC.plot <- ggplot(mod_comps, aes(x=desc, y=WAIC, color=scale)) +
                geom_point() +
                geom_line(aes(group=scale)) +
                theme_classic()
  ggsave(paste0(output_path,"WAIC_comparisons.pdf"),
         plot=WAIC.plot,width=80,height=90,units="cm")
  
  AUC_plot <- ggplot(mod_comps,aes(x=desc, y=AUC, color=scale)) +
              geom_point() +
              geom_line(aes(group=scale)) +
              theme_classic()


  ggsave(paste0(output_path,"Fit_comparisons.pdf"),
         plot=AUC_plot,width=80,height=90,units="cm")
  
  tiff(filename=paste0(output_path,"model_comp.tif"),
       width=18,height=9,units="cm",res=600)
  print(AUC_plot + WAIC.plot & theme(legend.position="none"))
  dev.off()
  