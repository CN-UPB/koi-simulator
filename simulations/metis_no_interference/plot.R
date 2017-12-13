#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

numMs <- c(20,30,40,50,60,70,80)

# Increase font size for all plots
#theme_update(base_size=24)

#delays_period_exploration(numMs,1,"./fixed_period_size_all_mcs/",9,
#                          figname="delays_fixed_period_fixed_size_all_mcs.pdf",
#                          tablename="delays_fixed_period_fixed_size_all_mcs.tab")
#deadlines_passed("deadlines_fixed_period_fixed_size_all_mcs.pdf",
#                 "delays_fixed_period_fixed_size_all_mcs.tab",
#                 "Deadline Overrunning for Burst Traffic with Fixed Packet Sizes")

#delays_period_exploration(numMs,1,"./fixed_period_size/",9,
#                          figname="delays_fixed_period_fixed_size",
#                          tablename="delays_fixed_period_fixed_size.tab")

deadlines_passed_fixed <- function(){
  figname<-"deadlines_fixed_period_fixed_size"
  tablename<-"delays_fixed_period_fixed_size.tab"
  figdir <- "./figures/"
  periods = c("1ms","2ms","5ms")
  data <- fread(paste0(figdir,tablename))
  # zoomPlot <- ggplotGrob(ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
  #   coord_cartesian(xlim = c(0.9, 1.1),ylim=c(0.9960,1.0),expand = FALSE)+
  #   geom_line(aes(group=period))+
  #   geom_point()+
  #   ylab("Pr")+
  #   xlab("#Devices")+
  #   theme_minimal()+
  #   theme(legend.position="none",plot.background = element_rect(colour = "black")))
  pl <- ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
#    annotation_custom(grob = zoomPlot, xmin = 3, xmax = 7, ymin = 0.62, ymax = 1.0)+
    geom_line(aes(group=period))+
    geom_point()+
    ylab("Pr[ packet delay < 1ms ]")+
    xlab("Nr. of Devices")+
    theme_bw(base_size=18)+
    theme(legend.justification=c(0,0), legend.position=c(0,0))
  jpeg(filename=paste0(figdir,figname,".jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(paste0(figdir,figname,".pdf"))
  plot(pl)
  dev.off()
}
deadlines_passed_fixed()

#delays_period_exploration(numMs,1,"./gaussian_period_fixed_size/",9,
#                          figname="delays_gaussian_period_fixed_size",
#                          tablename="delays_gaussian_period_fixed_size.tab")

deadlines_passed_gaussian <- function(){
  figname<-"deadlines_gaussian_period_fixed_size"
  tablename<-"delays_gaussian_period_fixed_size.tab"
  figdir <- "./figures/"
  periods = c("1ms","2ms","5ms")
  data <- fread(paste0(figdir,tablename))
  zoomPlot <- ggplotGrob(ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
                           coord_cartesian(xlim = c(1.8, 4.02),ylim=c(0.9925,1.001),expand = FALSE)+
                           geom_line(aes(group=period))+
                           geom_point()+
                           ylab("Pr")+
                           xlab("#Devices")+
                           theme_minimal()+
                           theme(legend.position="none",plot.background = element_rect(colour = "black")))
  pl <- ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
    geom_line(aes(group=period))+
    geom_point()+
    ylab("Pr[ packet delay < 1ms ]")+
    xlab("Nr. of Devices")+
    annotation_custom(grob = zoomPlot, xmin = 5, xmax = 7.4, ymin = 0.25, ymax = 0.9)+
    theme_bw(base_size=18)+
    theme(legend.justification=c(0,0), legend.position=c(0,0))
  jpeg(filename=paste0(figdir,figname),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,".pdf"))
  plot(pl)
  dev.off()
}
deadlines_passed_gaussian()

#delays_period_exploration(numMs,1,"./gaussian_period_and_size/",9,
#                          figname="delays_gaussian_period_and_size",
#                          tablename="delays_gaussian_period_and_size.tab")

deadlines_passed_gaussian <- function(){
  figname<-"deadlines_gaussian_period_and_size"
  tablename<-"delays_gaussian_period_and_size.tab"
  figdir <- "./figures/"
  periods = c("1ms","2ms","5ms")
  data <- fread(paste0(figdir,tablename))
  zoomPlot <- ggplotGrob(ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
                           coord_cartesian(xlim = c(1.8, 4.02),ylim=c(0.9925,1.001),expand = FALSE)+
                           geom_line(aes(group=period))+
                           geom_point()+
                           ylab("Pr")+
                           xlab("#Devices")+
                           theme_minimal()+
                           theme(legend.position="none",plot.background = element_rect(colour = "black")))
  pl <- ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
    geom_line(aes(group=period))+
    geom_point()+
    ylab("Pr[ packet delay < 1ms ]")+
    xlab("Nr. of Devices")+
    annotation_custom(grob = zoomPlot, xmin = 4, xmax = 7.2, ymin = 0.25, ymax = 0.95)+
    theme_bw(base_size=18)+
    theme(legend.justification=c(0,0), legend.position=c(0,0))
  jpeg(filename=paste0(figdir,figname),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,".pdf"))
  plot(pl)
  dev.off()
}
deadlines_passed_gaussian()