#### examine species predictions


now_area <- read.csv(paste0(output_path,"now_area_",short_modname,".csv"))
cc_area <- read.csv(paste0(output_path,"cc_area_",short_modname,".csv"))
stable_area <-  read.csv(paste0(output_path,"stable_area_",short_modname,".csv"))

# number of sp with <10 km area
bio.sp <- cc_area$sp[cc_area$bio<10]
length(bio.sp)

micro.sp <- cc_area$sp[cc_area$micro<10]
length(micro.sp)
micro.sp


# ABIEFRA preds
stable_area[stable_area$sp=="ABIEFRA",]

# PICERUB preds
stable_area[stable_area$sp=="PICERUB",]
