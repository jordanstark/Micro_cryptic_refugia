#### Setup ####
  # packages
    library(sp)
    library(raster)
    library(rasterVis)
    library(RColorBrewer)
    library(ggplot2)
    library(patchwork)
    
  fig_path <- paste0(output_path,"figures_",short_modname,"/")
  if(!dir.exists(fig_path)) dir.create(fig_path)

#### import data ####
  now.files <- list.files(paste0(model_inputs,"model_GIS/"),
                          full.names=T)
  cc.files <- list.files(paste0(model_inputs,"model_GIS_cc/"),
                         full.names=T)
  now.stack <- stack(now.files)
  names(now.stack) <- gsub("X_",names(now.stack),replacement="")
  cc.stack <- stack(cc.files)
  names(cc.stack) <- gsub("X_",names(cc.stack),replacement="")

#### extract layers for 'var' ####
  now <- now.stack[[grep("MAT",names(now.stack))]]
  cc <- cc.stack[[paste0(names(now),"_cc")]]
  
  if(!identical(paste0(names(now),"_cc"),names(cc))){
    stop("check raster order")
  }
  
  delta <- cc - now
  

  now.min <- min(cellStats(now,stat="min"))
  now.max <- max(cellStats(now,stat="max"))
  
  
  fine.bio.diff <- now[["fine_MAT"]] - now[["bio_MAT"]]
  micro.fine.diff <- now[["micro_MAT"]] - now[["fine_MAT"]]
  micro.bio.diff <- now[["micro_MAT"]] - now[["bio_MAT"]]
  
  micro.warming <- cc[["micro_MAT_cc"]] - now[["micro_MAT"]]
  
  max.diff <- max(cellStats(stack(fine.bio.diff,micro.fine.diff),"max"))
  min.diff <- min(cellStats(stack(fine.bio.diff,micro.fine.diff),"min"))
  
  ne.ext <- extent(290000,314000,3949000,3964000)
  
  fine.bio.section <- crop(fine.bio.diff,ne.ext)
  micro.fine.section <- crop(micro.fine.diff,ne.ext)
  
  
  scalebar_data <- data.frame(x=230000,xend=240000,y=3955000)
  scalebar_data_zoom <- data.frame(x=292000,xend=293000,y=3960500)

  
  
  pA <- gplot(now[["bio_MAT"]],maxpixels=1e7) +
          geom_tile(aes(fill=value)) +
          geom_segment(data=scalebar_data,aes(x=x,xend=xend,y=y,yend=y),
                       arrow=arrow(angle=90,ends="both",length=unit(0.01,"npc")))+
          geom_text(aes(x=mean(c(scalebar_data$x,scalebar_data$xend)),
                        y=scalebar_data$y+3000,label="10 km"),size=3) +
          theme_void() +
          scale_fill_distiller(name="MAT",
                               limits=c(now.min,now.max),
                               palette="RdBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("")  +
          coord_fixed(expand=F)
  
  pB <- gplot(now[["fine_MAT"]],maxpixels=1e7)+
          geom_tile(aes(fill=value)) +
          theme_void() +
          scale_fill_distiller(name="MAT",
                               limits=c(now.min,now.max),
                               palette="RdBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("")  +
          coord_fixed(expand=F)
  
  pC <- gplot(now[["micro_MAT"]],maxpixels=1e7)   +
          geom_tile(aes(fill=value)) +
          theme_void() +
          scale_fill_distiller(name="MAT",
                               limits=c(now.min,now.max),
                               palette="RdBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("")  +
          coord_fixed(expand=F)
  
  pD <- gplot(fine.bio.diff,maxpixels=1e7)+
          geom_tile(aes(fill=value)) +
          theme_void() +
          scale_fill_distiller(name="Change\nin MAT",
                               limits=c(min.diff,max.diff),
                               palette="RdYlBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("")  +
          coord_fixed(expand=F)
  

  pE <- gplot(micro.fine.diff,maxpixels=1e7)+
          geom_tile(aes(fill=value)) +
          theme_void() +
          scale_fill_distiller(name="Change\nin MAT",
                               limits=c(min.diff,max.diff),
                               palette="RdYlBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("") +
          coord_fixed(expand=F)

  pF <- gplot(micro.fine.section,maxpixels=1e7)+
          geom_tile(aes(fill=value)) +
          geom_segment(data=scalebar_data_zoom,aes(x=x,xend=xend,y=y,yend=y),
                       arrow=arrow(angle=90,ends="both",length=unit(0.01,"npc")))+
          geom_text(aes(x=mean(c(scalebar_data_zoom$x,scalebar_data_zoom$xend)),
                        y=scalebar_data_zoom$y+1000,label="1 km"),size=3) +
          theme_void() +
          scale_fill_distiller(name="Change\nin MAT",
                               limits=c(min.diff,max.diff),
                               palette="RdYlBu",
                               direction=-1,
                               na.value="white") +
          theme(legend.key.height=unit(0.5,"cm")) +
          xlab("") + ylab("")  +
          coord_fixed(expand=F) 
  
  
  sizes <- "
  AAAABBBBCCCCH
  AAAABBBBCCCCH
  AAAABBBBCCCCH
  DDDDEEEEFFFFH
  DDDDEEEEFFFFH
  DDDDEEEEFFFFH
  "

  
  tiff(filename=paste0(fig_path,"scale_clim.tif"),
       width=20,height=9,units="cm",res=600,compression="lzw")
  print(pA + pB + pC + pD + pE + pF + guide_area() +
          plot_layout(design=sizes,guides="collect") + 
          plot_annotation(tag_levels="A") &
          theme(plot.tag=element_text(size=10, hjust=-2,vjust=1.5),
                plot.tag.position=c(0,1)))
  dev.off()
  
  tiff(filename=paste0(fig_path,"micro_warming.tif"),
       width=18,height=9,units="cm",res=600,compression="lzw")
  print(gplot(micro.warming,maxpixels=1e7)+
          geom_tile(aes(fill=value)) +
          theme_void() +
          scale_fill_distiller(name="warming",
                               palette="RdYlBu",
                               direction=-1,
                               na.value="white") +
          xlab("") + ylab("")  +
          coord_fixed(expand=F) 
  )
  dev.off()
  
  
## comparison statistics
  cellStats(now,"mean")
  
  min.micro1km <- aggregate(now[["micro_MAT"]],fact=1000/30,fun=min,na.rm=T)
  max.micro1km <- aggregate(now[["micro_MAT"]],fact=1000/30,fun=max,na.rm=T)
  
  range.micro1km <- max.micro1km - min.micro1km

  fun.se <- function(x,na.rm=T) sd(x,na.rm=T)/sqrt(length(!is.na(x)))
  
  
  cellStats(range.micro1km,"mean")  
  cellStats(range.micro1km,"min") 
  cellStats(range.micro1km,"max") 
  
  
  min.fine1km <- aggregate(now[["fine_MAT"]],fact=1000/30,fun=min,na.rm=T)
  max.fine1km <- aggregate(now[["fine_MAT"]],fact=1000/30,fun=max,na.rm=T)
  
  range.fine1km <- max.fine1km - min.fine1km
  
  fun.se <- function(x,na.rm=T) sd(x,na.rm=T)/sqrt(length(!is.na(x)))
  
  
  cellStats(range.fine1km,"mean")  
  cellStats(range.fine1km,"min") 
  cellStats(range.fine1km,"max") 
  
  
  cellStats(micro.bio.diff,"max")
  cellStats(micro.bio.diff,"min")
  
  cellStats(micro.fine.diff,"max")
  cellStats(micro.fine.diff,"min")
  
  
  
  
  
  cellStats(now,"min")
  cellStats(now,"max")  
  
  cellStats(cc,"min")
  cellStats(cc,"max")  
  
  cellStats(micro.warming,"min")
  cellStats(micro.warming,"max")
  