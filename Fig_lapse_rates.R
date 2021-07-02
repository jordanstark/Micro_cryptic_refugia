#### Setup ####
  # packages
    library(lubridate)
    library(RColorBrewer)
    library(ggplot2)
    library(patchwork)
    
  fig_path <- paste0(output_path,"figures_",short_modname,"/")
  if(!dir.exists(fig_path)) dir.create(fig_path)

#### import data ####
  
  
  maxt_lapse <- read.csv(paste(climate_path,"maxtemp_synoptic_1900_2011.csv",sep=""))
  mint_lapse <- read.csv(paste(climate_path,"mintemp_synoptic_1900_2011.csv",sep=""))
  
  maxt_lapse$date <- mdy(maxt_lapse$date)
  mint_lapse$date <- mdy(mint_lapse$date)        
  maxt_lapse <- maxt_lapse[maxt_lapse$date < ymd("2000-01-01") & maxt_lapse$date > ymd("1969-12-31"),]
  mint_lapse <- mint_lapse[mint_lapse$date < ymd("2000-01-01") & mint_lapse$date > ymd("1969-12-31"),] 
  
  maxt_lapse$day <- yday(maxt_lapse$date)
  mint_lapse$day <- yday(mint_lapse$date)
  
  maxt_lapse$month <- month(maxt_lapse$date)
  mint_lapse$month <- month(mint_lapse$date)
  
  maxt_lapse$year <- year(maxt_lapse$date)
  mint_lapse$year <- year(mint_lapse$date)
  

#### model effect of intercept on lapse rate ####
  lm_maxt <- lm(slope ~ intercept, maxt_lapse)
  lm_mint <- lm(slope ~ intercept, mint_lapse)
  
  plot(slope ~ intercept,maxt_lapse,ylim=c(-0.022,0.021),xlim=c(-25,42))
  abline(a=lm_maxt$coefficients["(Intercept)"], b=lm_maxt$coefficients["intercept"], col="red")
  
  plot(slope ~ intercept,mint_lapse,ylim=c(-0.022,0.021),xlim=c(-25,42))
  abline(a=lm_mint$coefficients["(Intercept)"], b=lm_mint$coefficients["intercept"], col="red")
  
  maxt_slopeadj <- lm_maxt$coefficients["intercept"] * 4
  mint_slopeadj <- lm_mint$coefficients["intercept"] * 4
  
  
  
#### figure ####
  # lapse rate effect panels
  minx <- min(maxt_lapse$intercept,mint_lapse$intercept,na.rm=T)
  maxx <- max(maxt_lapse$intercept,mint_lapse$intercept,na.rm=T)
  miny <- min(maxt_lapse$slope,mint_lapse$slope,na.rm=T)
  maxy <- max(maxt_lapse$slope,mint_lapse$slope,na.rm=T)
  
  maxlapse_fig <- ggplot(maxt_lapse, aes(x=intercept,y=slope)) +
                    geom_point(shape=21,alpha=0.1) +
                    geom_smooth(method="lm",se=F,color="black") +
                    theme_classic() +
                    scale_x_continuous(limits=c(minx,maxx)) +
                    scale_y_continuous(limits=c(miny,maxy)) +
                    labs(title="maxima",
                         x="temperature at 0 m.a.s.l.",
                         y="change in temperature per m")
  minlapse_fig <- ggplot(mint_lapse, aes(x=intercept,y=slope)) +
                    geom_point(shape=21,alpha=0.1) +
                    geom_smooth(method="lm",se=F,color="black") +
                    theme_classic() +
                    scale_x_continuous(limits=c(minx,maxx)) +
                    scale_y_continuous(limits=c(miny,maxy),
                                       labels=NULL,name=NULL) +
                    labs(title="minima",
                         x="temperature at 0 m.a.s.l.")
    
  
  
  
  testdf <- data.frame(elev=rep(seq(300,2000,by=50),each=2))
  testdf$intercept <- rep(c(14,18))  
  testdf$temperature <- rep(c("cold","warm"))
  
  testdf$maxt_slope <- predict(lm_maxt,testdf)
  testdf$mint_slope <- predict(lm_mint,testdf)
  
  testdf$maxt <- testdf$maxt_slope*testdf$elev + testdf$intercept
  testdf$mint <- testdf$mint_slope*testdf$elev + testdf$intercept

  maxt_fig <- ggplot(testdf,aes(x=elev,y=maxt,color=temperature)) +
                geom_line() +
                theme_classic() + 
                scale_color_manual(values=c("blue","red")) +
                scale_y_continuous(limits=c(3,17))+
                geom_segment(aes(y=testdf$maxt[1]+0.5,
                                 yend=testdf$maxt[2]-0.5,
                                 x=300,
                                 xend=300),
                             color="black") +
                geom_segment(aes(y=testdf$maxt[69]+0.5,
                                 yend=testdf$maxt[70]-0.5,
                                 x=2000,
                                 xend=2000),
                             color="black") +
                geom_label(size=6/.pt,
                           aes(x=310,y=testdf$maxt[1]+0.75),
                          label=paste0(round(testdf$maxt[2]-testdf$maxt[1],1),
                                       "°C"),
                          label.size=NA,
                          color="black",
                          hjust=0) +
                geom_label(size=6/.pt,
                           aes(x=1950,y=testdf$maxt[70]-0.5),
                           label=paste0(round(testdf$maxt[70]-testdf$maxt[69],1),
                                        "°C"),
                           label.size=NA,
                           color="black",
                           hjust=1) +
                labs(x="elevation m",
                     y="temperature C") +
                theme(legend.position="none")
  
  mint_fig <- ggplot(testdf,aes(x=elev,y=mint,color=temperature)) +
                geom_line() +
                theme_classic() +
                scale_color_manual(values=c("blue","red")) +
                scale_y_continuous(limits=c(3,17),
                                   breaks=NULL,name=NULL) +
                geom_segment(aes(y=testdf$mint[1]+0.5,
                                 yend=testdf$mint[2]-0.5,
                                 x=300,
                                 xend=300),
                             color="black") +
                geom_segment(aes(y=testdf$mint[69]+0.5,
                                 yend=testdf$mint[70]-0.5,
                                 x=2000,
                                 xend=2000),
                             color="black") +
                geom_label(size=6/.pt,
                           aes(x=310,y=testdf$mint[1]+0.75),
                           label=paste0(round(testdf$mint[2]-testdf$mint[1],1),
                                        "°C"),
                           label.size=NA,
                           color="black",
                           hjust=0) +
                geom_label(size=6/.pt,
                           aes(x=1950,y=testdf$mint[70]-0.5),
                           label=paste0(round(testdf$mint[70]-testdf$mint[69],1),
                                        "°C"),
                           label.size=NA,
                           color="black",
                           hjust=1) +
                labs(x="elevation m",color="")
  
  
tiff(filename=paste0(fig_path,"lapse_effect.tif"),
     width=16.8,height=11,units="cm",res=600,compression="lzw")
  print(maxlapse_fig + minlapse_fig + maxt_fig + mint_fig +
    plot_layout(ncol=2) +
    plot_annotation(tag_levels="a") &
    theme(plot.tag=element_text(hjust=1,vjust=1.5),
          plot.tag.position=c(0,1),
          plot.title=element_text(hjust=0.5),
          text=element_text(size=6)))
dev.off()
