#### HMSC analysis on SVD ####

#### Setup ####
  # packages
    library(Hmsc)

    
  # paths
    # model_inputs <- "E:/Smokies_Veg/model_files/"
    # model_path <- "E:/Smokies_Veg/model_runs/"

#### import data ####
  presence_mat <- read.csv(paste0(model_inputs,"presence_data.csv"))
  presence_mat <- presence_mat[,3:length(presence_mat)]

  env_mat <- read.csv(paste0(model_inputs,"env_scaled.csv"))

#### specify and run models ####
  
  modlist <- vector(mode="list",length=length(formulas))
  
  names(modlist) <- names(formulas)

  thin <- 1000

  sampleFunc <- function(Hmsc.struct) {
    sampleMcmc(hM = Hmsc.struct,
               samples = 300,
               transient = 500*thin,
               thin = thin,
               nChains = 3,
               nParallel = 3)
  }
  
  
  for(i in 1:length(modlist)) {
    env_names_raw <- attr(terms.formula(formulas[[i]]),"term.labels")
    env_names <- sub(":",".",env_names_raw,fixed=T)
    
    #interaction may be in opposite order so we also need
    split_names <- strsplit(env_names,".",fixed=T)
    inv_names <- NA
    for(j in 1:length(split_names)){
      if(length(split_names[[j]])==1) {
        inv_names[j] <- split_names[[j]]
      } else if(length(split_names[[j]]==2)) {
        inv_names[j] <- paste0(split_names[[j]][2],".",split_names[[j]][1])
      } else stop("name error -- too many periods")
      
    }
    
    
    pa_envdat <- cbind(env_mat[,names(env_mat) %in% c(env_names,inv_names)],log_PlotSize=env_mat$log_PlotSize)

    if(length(pa_envdat) != length(env_names)+1) stop("error:some variables not found")
    
    pa_formula <- formula(paste0("~",paste0(names(pa_envdat),collapse=" + ")))
    
    pa.mod <- Hmsc(Y= presence_mat,
                    XData= pa_envdat,
                    XFormula= pa_formula,
                    XScale= F, #env vars already scaled
                    distr="probit")

    
    modlist[[i]] <- sampleFunc(pa.mod)

  }
  
  names(modlist) <- paste(names(formulas),rep("pa",length(names(formulas))),sep=".")
    
  save(modlist,file=paste0(model_path,outname))
  
  

 