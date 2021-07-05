# map distribution of single sp under different scenarios

library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

fig_sp_path <- paste0(output_path,"figures_",short_modname,"/sp_preds/")
if(!dir.exists(fig_sp_path)) dir.create(fig_sp_path)

fig_path <- paste0(output_path,"figures_",short_modname,"/")

parkbound <- readOGR(paste0(gis_path,"GRSM_BOUNDARY_POLYGON"))
parkbound <- parkbound[parkbound$OBJECTID<18,]


AUCs <- read.csv(paste0(output_path,"high_AUC_",short_modname,".csv"))


models <- data.frame(scale=rep(c("bio","fine","micro"),2),
                     scenario=rep(c("now","cc"),each=3))


### example northeast section for LIRITUL
ne.ext <- extent(290000,314000,3949000,3964000)

sp <- examplesp

bio.AUC <- round(AUCs$bio[AUCs$sp==sp],2)
fine.AUC <- round(AUCs$fine[AUCs$sp==sp],2)
micro.AUC <- round(AUCs$micro[AUCs$sp==sp],2)


rasterlist <- vector(mode="list",length=6)
names(rasterlist) <- paste(models$scale,models$scenario,sep=".")

for(j in 1:length(rasterlist)){
  scale <- models$scale[j]
  scen <- models$scenario[j]
  
  rasterlist[[j]] <- raster(paste0(output_path,"predicts_",short_modname,"/",
                                   scale,".pa",".",sp,".",scen,".tif"))
  
}

rstack <- stack(rasterlist,quick=T)
rstack <- crop(rstack,ne.ext)

bio.nodisp <- min(rstack$bio.now,rstack$bio.cc)
fine.nodisp <- min(rstack$fine.now,rstack$fine.cc)
micro.nodisp <- min(rstack$micro.now,rstack$micro.cc)

rstack <- stack(rstack,bio.nodisp,fine.nodisp,micro.nodisp)

area_km2 <- NA


for(j in 1:nlayers(rstack)){
  area <- cellStats(rstack[[j]],stat="sum")
  area_km2[j] <- round(area * 900 / 1000000)
  
}


plot_list <- list()

for(j in 1:nlayers(rstack)){
  pal <- ifelse(j<4,"BuGn",ifelse(j>6,"YlGn","PuBuGn"))
  
  plot_list[[j]] <- gplot(rstack[[j]],maxpixels=1e6) +
    geom_tile(aes(fill=value)) +
    theme_void() +
    scale_fill_distiller(name="",
                         limits=c(0,1),
                         palette=pal,
                         direction=1,
                         na.value="white") +
    theme(plot.title=element_text(hjust=0.5),
          text=element_text(size=10),
          axis.title.y=element_text(size=14),
          legend.key.height=unit(0.4,"cm")) +
    ylab("") +
    coord_fixed(expand=F) 
}

patch1 <- wrap_plots(plot_list[1:3],guides="collect") 

patch1[[1]] <- patch1[[1]] + 
  ggtitle(paste0("macroclimate")) + 
  ylab("historical")
patch1[[2]] <- patch1[[2]] + 
  ggtitle(paste0("interpolated"))
patch1[[3]] <- patch1[[3]] + 
  ggtitle(paste0("microclimate"))

patch2 <- wrap_plots(plot_list[4:6],guides="collect")
patch2[[1]] <- patch2[[1]] + ylab("dispersal")

patch3 <- wrap_plots(plot_list[7:9],guides="collect")
patch3[[1]] <- patch3[[1]] + ylab("stability")

full_patch <- patch1 / patch2 / patch3 


tiff(filename=paste0(fig_path,sp,"_NE.tif"),
     width=18,height=12,units="cm",res=600,compression="lzw")
print(full_patch)
dev.off()

### full figs for all sp

sp_table <- data.frame(sp=AUCs$sp,
                       AUC.bio=AUCs$bio,
                       AUC.fine=AUCs$fine,
                       AUC.micro=AUCs$micro,
                       stable.bio=NA,
                       stable.fine=NA,
                       stable.micro=NA,
                       disp.bio=NA,
                       disp.fine=NA,
                       disp.micro=NA)

for(i in 1:length(sp_table$sp)){ #
  sp <- sp_table$sp[i]
  bio.AUC <- round(sp_table$AUC.bio[i],2)
  fine.AUC <- round(sp_table$AUC.fine[i],2)
  micro.AUC <- round(sp_table$AUC.micro[i],2)
  fignum <- i + 7
  
  
  rasterlist <- vector(mode="list",length=6)
  names(rasterlist) <- paste(models$scale,models$scenario,sep=".")
  
  for(j in 1:length(rasterlist)){
    scale <- models$scale[j]
    scen <- models$scenario[j]
    
    rasterlist[[j]] <- raster(paste0(output_path,"predicts_",short_modname,"/",
                                     scale,".pa",".",sp,".",scen,".tif"))
    
  }
  
  rstack <- stack(rasterlist,quick=T)
  
  
  bio.nodisp <- min(rstack$bio.now,rstack$bio.cc)
  fine.nodisp <- min(rstack$fine.now,rstack$fine.cc)
  micro.nodisp <- min(rstack$micro.now,rstack$micro.cc)
  
  rstack <- stack(rstack,bio.nodisp,fine.nodisp,micro.nodisp)
  
  
  area_km2 <- NA
  
  
  for(j in 1:nlayers(rstack)){
    area <- cellStats(rstack[[j]],stat="sum")
    area_km2[j] <- round(area * 900 / 1000000)
    
  }
  
  sp_table$stable.bio[i] <- area_km2[7]/area_km2[1]
  sp_table$stable.fine[i] <- area_km2[8]/area_km2[2]
  sp_table$stable.micro[i] <- area_km2[9]/area_km2[3]
  
  sp_table$disp.bio[i] <- area_km2[4]/area_km2[1]
  sp_table$disp.fine[i] <- area_km2[5]/area_km2[2]
  sp_table$disp.micro[i] <- area_km2[6]/area_km2[3]
  
  plot_list <- list()
  
  for(j in 1:nlayers(rstack)){
    pal <- ifelse(j<4,"BuGn",ifelse(j>6,"YlGn","PuBuGn"))

    plot_list[[j]] <- gplot(rstack[[j]],maxpixels=1e6) +
      geom_tile(aes(fill=value)) +
      theme_void() +
      scale_fill_distiller(name="",
                           limits=c(0,1),
                           palette=pal,
                           direction=1,
                           na.value="white") +
      theme(plot.subtitle=element_text(hjust=0.5),
            plot.title=element_text(hjust=0.5)) +
      coord_fixed(expand=F) +
      labs(subtitle=paste0(area_km2[j]," km²"),y="")
  }

  patch1 <- wrap_plots(plot_list[1:3],guides="collect") 

  patch1[[1]] <- patch1[[1]] + 
    labs(title=paste0("macroclimate (AUC=",bio.AUC,")"),y="historical")
  patch1[[2]] <- patch1[[2]] + 
    labs(title=paste0("interpolated (AUC=",fine.AUC,")"))
  patch1[[3]] <- patch1[[3]] + 
    labs(title=paste0("microclimate (AUC=",fine.AUC,")"))
  
  patch2 <- wrap_plots(plot_list[4:6],guides="collect")
  patch2[[1]] <- patch2[[1]] + ylab("dispersal")
  
  patch3 <- wrap_plots(plot_list[7:9],guides="collect")
  patch3[[1]] <- patch3[[1]] + ylab("stability")
  
  full_patch <- patch1 / patch2 / patch3 + 
    plot_annotation(title=paste0("Fig S",fignum,". Probability of ",sp," occurrence")) &
    theme(axis.title.y=element_text(angle=90,size=14))
  

  tiff(filename=paste0(fig_sp_path,sp,".tif"),
       width=10,height=6,units="in",res=400,compression="lzw")
  print(full_patch)
  dev.off()
}


write.csv(sp_table,paste0(output_path,"sp_effects_table_",short_modname,".csv"))
