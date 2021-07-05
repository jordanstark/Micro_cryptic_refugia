#### Master script for climate scale SDM project ####
## Jordan Stark 
## Spring 2021

#########################################################
############# Setup and directory structure #############

## pathways
  scripts <- "E:/GithubRepos/Micro_cryptic_refugia/"
    # all scripts used in this master script  
  
  main <- "E:/Smokies_Veg/"
    # folder for each of the following paths
  setwd(main)
  
  db_path <- paste0(main,"csv_data/")
    # raw data exported from Smokies Veg Database
  
  gis_path <- paste0(main,"phys_GIS/")
    # park boundary and other physical GIS files
  
  fridley_GIS <- paste0(main,"phys_GIS/gsmnp_ascii/")
    # GIS files for calculating Fridley microclimate model
  
  climate_path <- paste0(main,"climate_data/")
    # raw and processed bioclim and microclimate rasters
    # and synoptic temperature estimates
  
  model_inputs <- paste0(main,"model_files/")
    # processed veg data for input to Hmsc
  
  model_path <- paste0(main,"model_runs/")
    # .RData files with saved Hmsc outputs
    # and scaling data
  
  output_path <- paste0(main,"outputs/")
    # model interpretation including prediction rasters
  
  
## required packages
  library(rstudioapi) # to run scripts in the background

  # all other packages are loaded in their own scripts
  ## GIS processing
    # raster (raster data)
    # sp (spatial points, also required for raster)
    # rgdal (spatial tranformation)
  ## climate calculation
    # ## dismo::biovars() was used as a template for calculating climate variables
      ## but it is not needed to run the scripts
  ## basic data processing
    # lubridate (for handling dates)
    # parallel (for parallel process)
    # tidyr (for pivoting long/wide data)
    # stringr (for text manipulation)
  ## Species modelling
    # Hmsc (Ovaskainen et al 2017)
  ## Data visualization
    # GGplot2 (for data figures)
    # rastervis (for GIS figures)
    # Patchwork (for multipanel figures)
  
  
## raw inputs and data sources
  ## macroclimate data
    # worldclim data at 0.5 degrees (downloaded if not in climate_path)
    # worldclim DEM in climate_path ("wc2.1_30s_elev")
  ## GRSM GIS data
    # "elev.txt" 30m park DEM in fridley_GIS -- from Fridley 2009
    # "GRSM_BOUNDARY_POLYGON" Park boundary in gis_path from GRSM datastore 
    # "tci.txt" 30m topographic convergence index in fridley_GIS
    # "totrad.txt" 30m annual radiation in fridley_GIS
    # "logsd.txt" 30m log-transformed stream distance in fridley_GIS
    # 365 daily solar radiation rasters (folder /rad/) in fridley_GIS
  ## microclimate model and lapse rates in climate_path
    # maxt model coefficients ("micro_mod/maxcoef_out.txt")
    # mint model coefficients ("micro_mod/mincoef_out.txt")
    # maxt lapse rates ("maxtemp_synoptic_1900_2011.csv") 
    # mint lapse rates ("mintemp_synoptic_1900_2011.csv") 
  ## species and plot data from Smokies Veg Database in db_path
    # "HerbdataSVD_3-14.csv" has all species records
    # "CarSpList_SVD.csv" has information about species ID
    # "SVD_w_spp_freq.csv" has growth havbit classifications
    # "Rasters_Mat21.csv" has spatial locations of all plots
  
##########################################################
########### Process climate data at all scales ###########
  
## check WorldClim estimates of downscaling and compare between GCMs
  jobRunScript(paste0(scripts,"Warming_exploration.R"),importEnv=T)
    # inputs: park boundary
    #         worldclim 2.1 data in climate_path
    #         worldclim downscaled GCM data in climate_path/GCM
    #         raster of elevations
    # outputs: supplemental figure comparing warming in each GCM
  
## download bioclim data, extract park, interpolate with DEM
  ## and add 4 degrees for climate change scenario
  
  jobRunScript(paste0(scripts,"Park_climate_wc21.R"),importEnv=T)
    # inputs: worldclim data at 0.5 degrees
    #         30m park DEM in fridley_GIS
    #         worldclim DEM in climate_path ("wc2.1_30s_elev")
    # outputs: 1km mean annual temp: "bio_MAT.tif" in climate_path
    #          1km max temp warmest month: "bio_MAXmo.tif" in climate_path
    #          1km min temp coldest month: "bio_MINmo.tif" in climate_path
    #          downscaled mean annual temp: "fine_MAT.tif" in climate_path
    #          downscaled max temp warmest month: "fine_MAXmo.tif" in climate_path
    #          downscaled min temp coldest month: "fine_MINmo.tif" in climate_path
    #          +4 degrees versions of all above files ("_plus4.tif") in climate_path
    #          rasters of easting and northing values ("Easting.tif", "Northing.tif") in gis_path
    
  
## calculate 1970-2000 microclimate and +4 degree microclimate
  ## and summarize to biovars (MAT, MAXmo, MINmo)
  jobRunScript(paste0(scripts,"Correct_TCI_raster.R"),importEnv=T)
    # inputs: "tci.txt" in fridley_GIS
    # outputs: "tci_cor.asc" in fridley_GIS (can be read as tci_cor.txt as well)
      # this has the '25' values (NA) replace with 95th quantile

  mainvars <- ls()
    #these are saved so that they aren't wiped out by the following scripts
  mainvars <- c(mainvars, "mainvars")
  
  source(paste0(scripts,"Micro_annual_invar.R"),local=T)
    # parallel fails with jobRunScript
    # so using source. Local=T means that the paths will copy
    # this takes ~ 15 mins to run
    # inputs: TCI ("tci.txt") in fridley_GIS
    #         totrad ("totrad.txt") in fridley_GIS
    #         elev ("elev.txt") in fridley_GIS
    #         strdist ("logsd.txt") in fridley_GIS
    #         365 daily solar radiation rasters (folder /rad/) in fridley_GIS
    #         maxt model coefficients ("micro_mod/maxcoef_out.txt") in climate_path
    #         mint model coefficients ("micro_mod/mincoef_out.txt") in climate_path
    # outputs: 365 rasters each in  climate_path '/tmpMin/' '/tmpMax/'
    # files are named as 'd_001.tif' ... 'd_365.tif' 
  

  source(paste0(scripts,"Micro_alldailyrasters.R"),local=T)
    # this takes ~ 6 hours
    # inputs: all of the above inputs from Micro_annual_invar.R 
    #         PLUS:
    #         maxt lapse rates ("maxtemp_synoptic_1900_2011.csv") in climate_path
    #         mint lapse rates ("mintemp_synoptic_1900_2011.csv") in climate_path
    #         365 daily rasters in climate_path generate as output of Micro_annualrasters.R
    # outputs: daily rasters with temperature records
    #         climate_path "/maxTs/" has 1970-2000 max temp
    #         climate_path "/minTs/" has 1970-2000 min temp
    #         climate_path "/maxTs_cc/" has +4 degees max temp
    #         climate_path "/minTs_cc/" has +4 degrees min temp
    #         named as 'y_1970_d_001.tif' ... 'y_1999_d_365.tif'
  
  source(paste0(scripts,"Micro_annual_summary.R"),local=T)
    # This takes ~1.5 hours
    # inputs: outputs of Micro_alldailyrasters.R in climate_path 
    #         "/maxTs/", "/minTs/", "/maxTs_cc/" and "/minTs_cc/"
    # outputs: rasters for each year for MAT, MAXmo and MINmo in
    #         climate_path "/summary/" for historical data
    #         climate_path "/summary_cc/" for +4 degrees scenario
  
  
  jobRunScript(paste0(scripts,"Micro_30y_summary.R"),importEnv=T)
    # this takes ~ 5 mins
    # inputs: outputs of "Micro_annual_summary.R" in climate_path
    #         "/summary/" for historical and "/summary_cc/" for +4 degrees
    # outputs: rasters of biovars 1,5,6 in climate_path
    #          historical: "micro_MAT.tif", "micro_MAXmo.tif", "micro_MINmo.tif"
    #          +4 degrees: "micro_MAT_plus4.tif", "micro_MAXmo_plus4.tif","micro_MINmo_plus4.tif"
  

  
##########################################################
######### Prep species and climate data for Hmsc #########
  
  jobRunScript(paste0(scripts,"Prep_points.R"),importEnv=T)
    # inputs: historical climate outputs of "Micro_30y_summary.R"
    #         fridley GIS files: "elev.txt","tci.txt","totrad.txt","logsd.txt"
    #         rasters of easting and northing from "Park_climate.R" 
    #         Vegetation database files:
    #           "HerbdataSVD_3-14.csv","CarSpList_SVD.csv","Rasters_May21.csv"
    # outputs: all files needed to run Hmsc models in model_inputs
    #         "presence_data.csv"
    #         "env_scaled.csv", "scalepars.csv"
  
##########################################################
######################## Run SDMs ########################
  # 'HMSC_models.R' takes outputs of "Prep_points.R" as inputs
  # and also requires specification of model formula and scales here
  # log(plotsize) is added automatically to the model formula
  # output is a list of length(scales) models
  # output is as a .RData file in "model_path"
  
  # these are parallelized internally 3x - so can run 5 together
  # each model takes about 20h, so 2.5 days for all three
  

  outname <- "MATxtopo+space.RData"
  formulas <- list(micro= ~ (micro_MAT + micro_MAT2)*(log_tci + totrad)+Easting+Northing,
                   fine=  ~ (fine_MAT + fine_MAT2)*(log_tci + totrad)+Easting+Northing,
                   bio=   ~ (bio_MAT + bio_MAT2)*(log_tci + totrad)+Easting+Northing)
  jobRunScript(paste0(scripts,"HMSC_models.R"),importEnv=T)
  
  outname <- "MAT+topo+space.RData"
  formulas <- list(micro= ~ (micro_MAT + micro_MAT2)+(log_tci + totrad)+Easting+Northing,
                   fine=  ~ (fine_MAT + fine_MAT2)+(log_tci + totrad)+Easting+Northing,
                   bio=   ~ (bio_MAT + bio_MAT2)+(log_tci + totrad)+Easting+Northing)
  jobRunScript(paste0(scripts,"HMSC_models.R"),importEnv=T)
  
  outname <- "elev+space.RData"
  formulas <- list(elev= ~ elev + Easting + Northing)
  jobRunScript(paste0(scripts,"HMSC_models.R"),importEnv=T)
  

##########################################################
################## model goodness of fit #################
## overall model comparison
  
  # list of models to compare (shortened name = .RData file)
  models <- list("elev"="elev+space.RData",
                 "MAT+topo"="MAT+topo+space.RData",
                 "MATxtopo"="MATxtopo+space.RData"
                 )
  jobRunScript(paste0(scripts,"HMSC_model_comparison.R"),importEnv=T)
    # inputs: .RData files specified above in 'models'
    #         each .RData file is created by one run of"HMSC_models.R"
    #         and contains presence/absence
    #         usually for the three sets of climate data (for elev only, just one model)
    # outputs: "model_comparisons.csv" in output_path has 
    #                 WAIC, AUC, psrf for each model
    #         "WAIC_comparisons.pdf" in output_path 
    #         "Fit_comparisons.pdf" in output_path
  
## compare AUC by species within a set of models (3 min each, can run together)


    modlist_path <- "MATxtopo+space.RData"
    short_modname <- "MATxtopo"
    jobRunScript(paste0(scripts,"HMSC_species_fit.R"),importEnv=T)
    # inputs: models in modlist_path
    # outputs: list of sp with AUC >0.8 in all mods
    #          "high_AUC"_short_modname in output_path
    
    modlist_path <- "MAT+topo+space.RData"
    short_modname <- "MAT+topo"
    jobRunScript(paste0(scripts,"HMSC_species_fit.R"),importEnv=T)
    


##########################################################
############# Predict presence and abundance #############
## prep rasters that are used to make predictions
  jobRunScript(paste0(scripts,"HMSC_predict_prep_GIS.R"),importEnv=T)
  # inputs: historical and cc climate outputs of "Micro_30y_summary.R"
  #         fridley GIS files: "elev.txt","tci.txt","totrad.txt","logsd.txt"
  #         rasters of easting and northing from "Park_climate.R" 
  # outputs: raster stacks with all raw files for Hmsc models in model_files
  #          "model_GIS/" has historical, and "model_GIS_cc/" has +4 degrees
  
  mainvars <- ls()
  #these are saved so that they aren't wiped out by the following scripts
  mainvars <- c(mainvars, "mainvars")
  
  # creating the prediction rasters takes ~30 mins per model

  modlist_path <- "MATxtopo+space.RData"
  short_modname <- "MATxtopo"
  source(paste0(scripts,"HMSC_species_predict.R"))
  # inputs: rasters from "HMSC_predict_prep_GIS.R"
  #         model list in modlist_path generated by "HMSC_models.R"
  #         high AUC species from "HMSC_species_fit.R" in output_path
  # outputs: rasters of probability of occurence
  #          for scenarios 'now' and 'cc'
  #          in folder output_path/predicts_short_modname


  ## prediction rasters for MATxtopo+space model for all species, regardless of AUC
  modlist_path <- "MATxtopo+space.RData"
  short_modname <- "MATxtopo"
  source(paste0(scripts,"HMSC_all_species_predict.R"))
  
  
  
##########################################################
################ Compare species habitat #################
  #~ 15-30 mins, can run at same time
  
  short_modname <- "MATxtopo"
  jobRunScript(paste0(scripts,"Summarize_distributions.R"),importEnv=T)
  # inputs: rasters from 'HMSC_species_predict.R' in predicts_folder
  # outputs: number of sp in each scenario in outputs_folder/spsummary_short_modname
  #          fraction stable sp at each scale
  #          csv files with area occupied in each scenario and stable area
  #          rasters of stable area in outputs_folder / stable_short_modname
  


  
  short_modname <- "MATxtopo_all"
  jobRunScript(paste0(scripts,"Summarize_distributions.R"),importEnv=T)
  
##########################################################
############ Make figures comparing results ##############
  modlist_path <- "MATxtopo+space.RData"
  short_modname <- "MATxtopo"

# Climate figures
  jobRunScript(paste0(scripts,"Fig_micro_minus_macro.R"),importEnv=T)

# sample points, species change and hex figures
  jobRunScript(paste0(scripts,"fig_sp_topo.R"),importEnv=T)
  
# lapse rate change with climate figure
  jobRunScript(paste0(scripts,"Fig_lapse_rates.R"),importEnv=T)

# parkwide change in species distribution 
  # (1 tif per sp of whole park, plus zoomed in fig for example sp)
  
  examplesp <- "LIRITUL"
  jobRunScript(paste0(scripts,"Fig_sp_abundance.R"),importEnv=T)
  

# check species outputs
  # Examine_sp_outputs.R
  

# repeat figures for all species
  modlist_path <- "MATxtopo+space.RData"
  short_modname <- "MATxtopo_all"
  
  