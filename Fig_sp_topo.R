#### Setup ####
  # packages
    library(sp)
    library(raster)
    library(rgdal)
    library(ggplot2)
    library(RColorBrewer)
    library(rasterVis)
    library(tidyr)    
    library(patchwork)
    library(ggrepel)
    library(Hmsc)

fig_path <- paste0(output_path,"figures_",short_modname,"/")
if(!dir.exists(fig_path)) dir.create(fig_path)


#### import data ####
  # raster layers
    now.files <- list.files(paste0(model_inputs,"model_GIS/"),
                            full.names=T)
    cc.files <- list.files(paste0(model_inputs,"model_GIS_cc/"),
                           full.names=T)
    now.stack <- stack(now.files)
    names(now.stack) <- gsub("X_",names(now.stack),replacement="")
    cc.stack <- stack(cc.files)
    names(cc.stack) <- gsub("X_",names(cc.stack),replacement="")
    
    log_strdist <- raster(paste0(fridley_GIS,"logsd.txt"))
    crs(log_strdist) <- CRS("+proj=utm +zone=17 +datum=NAD27")
    log_strdist <- crop(log_strdist,now.stack)
      
  # models
    frac.bio.stable <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                     "frac.bio.stable.tif"))
    frac.fine.stable <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                     "frac.fine.stable.tif"))
    frac.micro.stable <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                     "frac.micro.stable.tif"))
    
    
    num.bio.now <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                    "num.bio.now.tif"))
    num.fine.now <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                 "num.fine.now.tif"))
    num.micro.now <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                 "num.micro.now.tif"))
   
    num.now <- stack(num.bio.now,num.fine.now,num.micro.now)
    
    num.bio.cc <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                 "num.bio.cc.tif"))
    num.fine.cc <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                  "num.fine.cc.tif"))
    num.micro.cc <- raster(paste0(output_path,"spsummary_",short_modname,"/",
                                   "num.micro.cc.tif"))
    
    num.cc <- stack(num.bio.cc,num.fine.cc,num.micro.cc)
    
    
    change.cc <- num.cc/num.now
    names(change.cc) <- c("frac.bio.expand","frac.fine.expand","frac.micro.expand")
    
    
  # park boundary to sample
    parkbound <- readOGR(paste(gis_path,"GRSM_BOUNDARY_POLYGON",sep=""))
    parkbound <- spTransform(parkbound,crs(now.stack))
    parkbound <- parkbound[parkbound$OBJECTID<18,]
    
    
  # calc AUC values
    load(paste0(model_path,modlist_path))
    
    
    pa.mods <- modlist[grep("pa",names(modlist))]
    
    calcAUC <- function(mod) {
      pred <- computePredictedValues(mod)
      AUC <- evaluateModelFit(mod,pred)[["AUC"]]
    }
    
    AUC <- data.frame(sapply(pa.mods,calcAUC))
    AUC <- data.frame(AUC)
    
    AUCdf <- cbind(sp=modlist[[1]]$spNames,AUC)
    names(AUCdf)[2:length(AUCdf)] <- simplify2array(strsplit(names(AUCdf[2:length(AUCdf)]),".",fixed=T))[1,]
    

#### extract random points ####
    set.seed(3284)
    pts <- spsample(parkbound,10000,type="regular")
    
    allstack <- stack(now.stack,
                      cc.stack[[grep("cc",names(cc.stack))]],
                      log_strdist,
                      frac.bio.stable,
                      frac.fine.stable,
                      frac.micro.stable,
                      change.cc)
    
    ptdat <- data.frame(raster::extract(allstack,pts))
    ptdat$ptID <- 1:length(ptdat[,1])

 
    
    ptdat$micro_MAT_delta <- ptdat$micro_MAT_cc - ptdat$micro_MAT
    ptdat$fine_MAT_delta <- ptdat$fine_MAT_cc - ptdat$fine_MAT
    ptdat$bio_MAT_delta <- ptdat$bio_MAT_cc - ptdat$bio_MAT
    
    
    
    plotlist <- vector(mode="list",length=12)
    
    scales <- c("bio","fine","micro")
    scale_names <- c("macroclimate","interpolation","microclimate")
    
    nbins <- 10
    
    sdlabs <- c(10,100,1000)
    
    sdbreaks <- log(sdlabs)
    
    
    mean10 <- function(x) ifelse(sum(!is.na(x))>=10, mean(x,na.rm=T),NA)
    xlim <- c(2.2,7)
    ylim <- c(250,2000)

    
    p1 <- ggplot(ptdat,aes(x=logsd,y=elev,z=bio_MAT)) +
                  ylim(ylim) +
                  ggtitle("macroclimate") 
    p2 <- ggplot(ptdat,aes(x=logsd,y=elev,z=fine_MAT)) +
                  scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL) +
                  ggtitle("interpolation")
    p3 <- ggplot(ptdat,aes(x=logsd,y=elev,z=micro_MAT)) +
                  scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL) +
                  ggtitle("microclimate")
      
      
    row1 <-  (p1|p2|p3) + 
      plot_layout(guides="collect")
      
    
    row1 <- row1 &
              stat_summary_hex(aes(color=stat(value)),
                               fun="mean10",bins=nbins) &
              scale_fill_distiller(aesthetics=c("color","fill"),
                                   name="MAT        ",
                                   limits=c(6,17),
                                   palette="RdBu",direction=-1) &
              scale_x_continuous(limits=xlim,labels=NULL,breaks=NULL)

    p4 <- ggplot(ptdat,aes(x=logsd,y=elev,z=bio_MAT_delta)) +
                  ylim(ylim)
    p5 <- ggplot(ptdat,aes(x=logsd,y=elev,z=fine_MAT_delta)) +
                  scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)
    p6 <- ggplot(ptdat,aes(x=logsd,y=elev,z=micro_MAT_delta)) +
                  scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)
    
    row2 <- (p4 | p5 | p6) + 
      plot_layout(guides="collect")
    
    row2 <- row2 &
            stat_summary_hex(aes(color=stat(value)),
                             fun="mean10",bins=nbins) &
            scale_fill_distiller(aesthetics=c("color","fill"),
                                 name="warming ",
                                 limits=c(2.1,4.1),
                                 palette="RdYlBu",direction=-1) &
            scale_x_continuous(limits=xlim,labels=NULL,breaks=NULL)
    

    
    p7 <- ggplot(ptdat, aes(x=logsd, y=elev, z=frac.bio.stable)) +
            ylim(ylim)
    p8 <- ggplot(ptdat, aes(x=logsd, y=elev, z=frac.fine.stable))+
            scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)
    p9 <- ggplot(ptdat, aes(x=logsd, y=elev, z=frac.micro.stable))+
            scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)

    
    row3 <- (p7|p8|p9) + 
      plot_layout(guides="collect")
    
    row3 <- row3 &                    
                    stat_summary_hex(aes(color=stat(value)),
                                     fun="mean10",bins=nbins) &
                    scale_fill_distiller(aesthetics=c("color","fill"),
                                         name="stability   ",
                                         limits=c(0.13,0.7),
                                         palette="YlGn",direction=1) &

                    scale_x_continuous(limits=xlim, labels=NULL,breaks=NULL)

    
    p10 <- ggplot(ptdat, aes(x=logsd,y=elev,z=frac.bio.expand)) +
            ylim(ylim)
    p11 <- ggplot(ptdat, aes(x=logsd,y=elev,z=frac.fine.expand)) +
            scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)
    p12 <- ggplot(ptdat, aes(x=logsd,y=elev,z=frac.micro.expand)) +
            scale_y_continuous(limits=ylim,labels=NULL,breaks=NULL)
    
    
    row4 <- (p10|p11|p12) + 
      plot_layout(guides="collect")
    
    row4 <- row4 &  
                  stat_summary_hex(aes(color=stat(value)),fun="mean10",bins=nbins) &
                  scale_fill_distiller(aesthetics=c("color","fill"),
                                       name="dispersal",
                                       limits=c(0.7,2.7),
                                       palette="PuBuGn",direction=1) &
                  scale_x_continuous(limits=xlim,breaks=sdbreaks,labels=sdlabs)
    
    
    
    
    full.patch <- row1 / row2 / row3 / row4 
    
    
    
    full.patch <- full.patch &
                    theme_classic() &
                    xlab("") &
                    ylab("") &
                    theme(plot.title=element_text(hjust=0.5,face="bold"))
    
    layout <- "
                #AAAAAAAAAAAAAAAA
                #AAAAAAAAAAAAAAAA
                #AAAAAAAAAAAAAAAA
                #AAAAAAAAAAAAAAAA
                #BBBBBBBBBBBBBBBB
                EBBBBBBBBBBBBBBBB
                EBBBBBBBBBBBBBBBB
                EBBBBBBBBBBBBBBBB
                ECCCCCCCCCCCCCCCC
                ECCCCCCCCCCCCCCCC
                ECCCCCCCCCCCCCCCC
                #CCCCCCCCCCCCCCCC
                #DDDDDDDDDDDDDDDD
                #DDDDDDDDDDDDDDDD
                #DDDDDDDDDDDDDDDD
                #DDDDDDDDDDDDDDDD
                ######FFFFFF#####"
    
    final.fig <- full.patch + 
                  grid::textGrob("elevation (m)",rot=90,vjust=1,
                                 gp=grid::gpar(fontface="bold")) + 
                  grid::textGrob("stream distance (m)", vjust=-1,
                                 gp=grid::gpar(fontface="bold")) + 
                  plot_layout(design=layout)


    #final.fig
    
    tiff(filename=paste0(fig_path,"scale_topo_warming.tif"),
         width=18,height=20,units="cm",res=600)
    print(final.fig)
    dev.off()
 
    
    
# map figure

    diff.stable.mf <- frac.micro.stable - frac.fine.stable
    diff.stable.fb <- frac.fine.stable - frac.bio.stable
    diff.stable.mb <- frac.micro.stable - frac.bio.stable
    
    diff.disp.mf <- change.cc$frac.micro.expand - change.cc$frac.fine.expand
    diff.disp.fb <- change.cc$frac.fine.expand - change.cc$frac.bio.expand
    diff.disp.mb <- change.cc$frac.micro.expand - change.cc$frac.bio.expand
    
    scale_cols <- c("tomato","goldenrod1","royalblue3")
    
    
    stable.mf <- gplot(diff.stable.mf,maxpixels=1e7)+
                  geom_tile(aes(fill=value)) +
                  theme_void() +
                  scale_fill_gradient2(name="difference\nin stability",
                                       low=scale_cols[2],
                                       high=scale_cols[3],
                                       na.value="white")  +
                  theme(text=element_text(size=12),
                        plot.title=element_text(hjust=0.5)) +
                  coord_fixed(expand=F)

    
    tiff(filename=paste0(fig_path,"cryptic_refugia.tif"),
         width=18,height=9,units="cm",res=600)
    print(stable.mf)
    dev.off()
    
    ne.ext <- extent(290000,314000,3949000,3964000)
    
    section <- crop(diff.stable.mf,ne.ext)
    
    stable.mf.section <-gplot(section,maxpixels=1e7)+
                              geom_tile(aes(fill=value)) +
                              theme_void() +
                              scale_fill_gradient2(name="difference\nin stability",
                                                   low=scale_cols[2],
                                                   high=scale_cols[3],
                                                   na.value="white",
                                                   limits=c(cellStats(diff.stable.mf,"min"),
                                                            cellStats(diff.stable.mf,"max"))) +
                              theme(text=element_text(size=12),
                                    plot.title=element_text(hjust=0.5)) +
                              coord_fixed(expand=F)
    
    
    tiff(filename=paste0(fig_path,"cryptic_refugia_ne.tif"),
         width=18,height=9,units="cm",res=600)
    print(stable.mf.section)
    dev.off()
    
    
    des <- "AAAABBC"
    
    cryptic_fig <- stable.mf + stable.mf.section +
                    plot_layout(guides="collect") +
                    plot_annotation(tag_levels="A")
    
    tiff(filename=paste0(fig_path,"cryptic_refugia_fig.tif"),
         width=18,height=9,units="cm",res=600)
    print(cryptic_fig)
    dev.off()
    
    
    stable.mb <- gplot(diff.stable.mb,maxpixels=1e7)+
                  geom_tile(aes(fill=value)) +
                  theme_void() +
                  scale_fill_gradient2(name="difference\nin stability",
                                       low=scale_cols[1],
                                       high=scale_cols[3],
                                       na.value="white")  +
                  theme(text=element_text(size=12),
                        plot.title=element_text(hjust=0.5)) +
                  coord_fixed(expand=F)
    
    
    mb.section <- crop(diff.stable.mb,ne.ext)
    
    stable.mb.section <-gplot(mb.section,maxpixels=1e7)+
                          geom_tile(aes(fill=value)) +
                          theme_void() +
                          scale_fill_gradient2(name="difference\nin stability",
                                               low=scale_cols[1],
                                               high=scale_cols[3],
                                               na.value="white",
                                               limits=c(cellStats(diff.stable.mb,"min"),
                                                        cellStats(diff.stable.mb,"max"))) +
                          theme(text=element_text(size=12),
                                plot.title=element_text(hjust=0.5)) +
                          coord_fixed(expand=F)
    
    tiff(filename=paste0(fig_path,"micro_refugia.tif"),
         width=18,height=9,units="cm",res=600)
    print(stable.mb)
    dev.off()
    
    tiff(filename=paste0(fig_path,"micro_refugia_ne.tif"),
         width=18,height=9,units="cm",res=600)
    print(stable.mb.section)
    dev.off()

    
    
    stable_area <- read.csv(paste0(output_path,"stable_area_frac",short_modname,".csv"))
    change_area <- read.csv(paste0(output_path,"change_area_",short_modname,".csv"))
    now_area <- read.csv(paste0(output_path,"now_area_",short_modname,".csv"))
    cc_area <- read.csv(paste0(output_path,"cc_area_",short_modname,".csv"))
    
    identical(now_area$sp,cc_area$sp)
    
    ratio_change <- cc_area[,3:5] / now_area[,3:5]
    ratio_change <- cbind(sp=now_area$sp,ratio_change)
    
    stable_cleaned <- stable_area[stable_area$micro<0.98 & 
                                    stable_area$fine<0.98 &
                                    stable_area$bio<0.98,]
    
    
    stable_long <- pivot_longer(stable_area,
                                cols=c("fine","bio","micro"),
                                names_to="modtype",
                                values_to="frac_area_stable")
    stable_cleaned_long <-  pivot_longer(stable_cleaned,
                                         cols=c("fine","bio","micro"),
                                         names_to="modtype",
                                         values_to="frac_area_stable")
    stable_cleaned_long <- stable_cleaned_long[stable_cleaned_long$sp != "ATHYFIL",]
    
    

    

    change_long <- pivot_longer(change_area,
                                cols=c("fine","bio","micro"),
                                names_to="modtype",
                                values_to="change_area_km2")
    
    change_cleaned_long <- change_long[change_long$sp %in% stable_cleaned_long$sp,]
    
    
    ratio_change_long <- pivot_longer(ratio_change,
                                      cols=c("fine","bio","micro"),
                                      names_to="modtype",
                                      values_to="future_current")
    ratio_cleaned_long <- ratio_change_long[ratio_change_long$sp  %in% stable_cleaned_long$sp,]
    
    
    stable_plot <- ggplot(stable_cleaned_long,aes(x=modtype, y=frac_area_stable)) +
                      geom_violin(aes(fill=modtype,group=modtype),show.legend=F,alpha=0.5) +
                      geom_point(size=0.4) +
                      geom_line(aes(group=sp),alpha=0.6,size=0.2) +
                      geom_hline(yintercept=1,color="grey30",
                                 size=0.8,linetype="dashed") +
                      scale_fill_discrete(type=scale_cols) +
                      theme_classic() +
                      xlab("") +
                      ylab("future area:current area\n") +
                      theme(text=element_text(size=10,color="black"),
                            axis.line = element_line(color="black",size=0.5),
                            axis.ticks = element_line(color="black",size=0.5)) +
                      scale_x_discrete(labels=c("macroclimate",
                                                "interpolation",
                                                "microclimate")) +
                      scale_y_continuous(expand=c(0,0),limits=c(0,1.05)) +
                      geom_text_repel(aes(label=sp),max.overlaps=5,
                                      min.segment.length=10,
                                      size=1.8,fontface="bold") 
    
    # change_plot <- ggplot(change_cleaned_long,aes(x=modtype,y=change_area_km2)) +
    #   geom_violin(aes(fill=modtype,group=modtype),show.legend=F) +
    #   geom_point() +
    #   geom_line(aes(group=sp),alpha=0.5) +
    #   scale_fill_discrete(type=scale_cols) +
    #   theme_classic() +
    #   xlab("") +
    #   ylab("change in area (km2)\n") +
    #   theme(text=element_text(size=12),
    #         axis.text=element_text(color="black")) +
    #   scale_x_discrete(labels=c("macroclimate","interpolation","microclimate")) +
    #   geom_text_repel(aes(label=sp),max.overlaps=4,min.segment.length=10,
    #                   size=2,fontface="bold") 
    
    
    ratio_plot <- ggplot(ratio_cleaned_long,aes(x=modtype,y=future_current)) +
                    geom_violin(aes(fill=modtype,group=modtype),show.legend=F,alpha=0.5) +
                    geom_point(size=0.4) +
                    geom_line(aes(group=sp),alpha=0.6,size=0.2) +
                    geom_hline(yintercept=1,color="grey30",
                               size=0.8,linetype="dashed") +
                    scale_fill_discrete(type=scale_cols) +
                    theme_classic() +
                    xlab("") +
                    ylab("") +
                    theme(text=element_text(size=10,color="black"),
                          axis.line = element_line(color="black",size=0.5),
                          axis.ticks = element_line(color="black",size=0.5)) +
                    scale_x_discrete(labels=c("macroclimate",
                                              "interpolation",
                                              "microclimate")) +
                    scale_y_continuous(expand=c(0,0)) +
                    geom_text_repel(aes(label=sp),max.overlaps=5,
                                    min.segment.length=10,
                                    size=1.8,fontface="bold") 
    
    tiff(filename=paste0(fig_path,"species_response.tif"),
         width=18,height=9,units="cm",res=600)
    print(stable_plot + ratio_plot + plot_annotation(tag_levels="A"))
    dev.off()
    
    
    
    spdat <- read.csv(paste(db_path,"HerbdataSVD_3-14.csv",sep=""))
    plotdat <- read.csv(paste(db_path,"Rasters_May21.csv",sep=""))[,c("PLOT_ID","DEM10")]


      
    sp_types <- read.csv(paste(db_path,"SVD_w_spp_freq.csv",sep=""),na.strings="")
    sp_types <- sp_types[,c("GEN.SPP","Growth.Habit")]
    sp_types <- sp_types[complete.cases(sp_types),]
    sp_types <- merge(sp_types,spdat[,c("NC_CODE","GEN.SPP")])
    sp_types$GEN.SPP <- NULL
    names(sp_types) <- c("Habit","sp")
    sp_types <- unique(sp_types)

    sp_types$TreeShrub <- F
    sp_types$TreeShrub[grep("Tree|Shrub",sp_types$Habit)] <- T
    herbtypes <- c("Forb|Herb|Graminoid|Vine")
    sp_types$Herb <- F
    sp_types$Herb[grep(herbtypes,sp_types$Habit)] <- T


    sp_types$SimpleHab <- ifelse(sp_types$TreeShrub==T & sp_types$Herb==F, "TreeShrub",
                                   ifelse(sp_types$Herb==T,"Herb",NA))

    alldb <- merge(spdat,plotdat)

    spmin <- aggregate(DEM10~NC_CODE,
                             alldb,
                             FUN=min,na.rm=T)

    names(spmin) <- c("sp","min.elev")

    spmean <- aggregate(DEM10~NC_CODE,
                         alldb,
                         FUN=mean,na.rm=T)
      
    names(spmean) <- c("sp","mean.elev")
      
    spvals <- merge(spmin,spmean)
    spsummary2 <- merge(spvals,sp_types,all.x=T)

    stable_long2 <- merge(stable_long,spsummary2)

      
    ggplot(stable_long2, aes(x=min.elev,y=frac_area_stable,color=modtype)) +
        geom_point() +
        theme_classic()
      
    tiff(filename=paste0(fig_path,"elev_stability.tif"),
         width=18,height=9,units="cm",res=600)
    print(ggplot(stable_long2, aes(x=mean.elev,y=frac_area_stable,color=modtype,color=modtype)) +
        geom_point(size=1.5,alpha=0.8) +
        geom_point(size=1.5,shape=1,stroke=0.5,alpha=0.8) +
        scale_color_discrete(type=scale_cols,
                             labels=c("macroclimate","interpolation","microclimate")) +
        geom_smooth(se=F,alpha=0.8) +
        theme_classic() +
        xlab("\nmean elevation (m)") +
        ylab("stable future:current area\n")  +    
        theme(text=element_text(size=10),
              axis.text=element_text(color="black"),
              legend.position=c(0.85,0.8),) +
        ylim(0,1) +
        labs(color=""))
    dev.off()
      
      AUC_long <- pivot_longer(AUCdf,micro:bio,
                               names_to="modtype",
                               values_to="AUC")
      stable_long3 <- merge(spsummary2, AUC_long)
      
      
      ggplot(stable_long3, aes(x=mean.elev,y=AUC,color=modtype)) +
        geom_point(size=2) +
        scale_color_discrete(type=scale_cols,
                             labels=c("macroclimate","interpolation","microclimate")) +
        theme_classic() +
        xlab("\nmean elevation (m)") +
        ylab("AUC\n")  +    
        theme(text=element_text(size=14),
              axis.text=element_text(color="black"),
              legend.text=element_text(size=12)) +
        labs(color="scale")
  
  tiff(filename=paste0(fig_path,"elev_AUC.tif"),
       width=18,height=9,units="cm",res=600)
      print(ggplot(stable_long3, aes(x=mean.elev,y=AUC,color=modtype,shape=SimpleHab)) +
              geom_point(size=2) +
              geom_hline(yintercept=0.8,color="grey30",size=0.8,linetype="dashed") +
              scale_color_discrete(type=scale_cols,aesthetics=c("color","fill"),
                                   labels=c("macroclimate","interpolation","microclimate")) +
              theme_classic() +
              xlab("\nmean elevation (m)") +
              ylab("AUC\n")  +    
              theme(text=element_text(size=10,color="black"),
                    axis.line = element_line(color="black",size=0.5),
                    axis.ticks = element_line(color="black",size=0.5),
                    legend.position=c(0.85,0.2)) +
              guides(color="legend",shape="none") +
              labs(color="",shape="Habit"))
    dev.off()
      
    