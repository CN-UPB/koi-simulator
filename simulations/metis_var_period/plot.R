#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

numMs <- c(20,30,40,50,60,70,80)
numCells <- 9

#delays_period_exploration(numMs,10,"./full_results/",numCells)

deadlines_passed <- function(figname,tablename){
  figdir <- "./figures/"
  jpeg(filename=paste0(figdir,figname),quality=100,type="cairo")
  periods = c("1ms","2ms","5ms")
  data <- fread(paste0(figdir,tablename))
  pl <- ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
    geom_line(aes(group=period))+
    geom_point()+
    ylab("Pr[ packet delay < 1ms ]")+
    xlab("Nr. of Devices")+
    theme_bw(base_size=18)+
    theme(legend.justification=c(0,0), legend.position=c(0,0))
  jpeg(filename=paste0(figdir,figname,".jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,".pdf"))
  plot(pl)
  dev.off()
}
deadlines_passed("deadlines_full_spectrum",
                 "delays_smaller_1ms_by_ms.tab")

# 10 ms Runs
# resDir <- "./results_10ms/"
# figname <- "full_spectrum"
# tablename <- "full_spectrum_10ms.tab"
# figdir <- paste0(resDir,"../figures/")
# periods = c("2ms","5ms","10ms")
# numMS <- numMs
# ecdfs <- data.table()
# numConfs <- length(numMS)
# runs <- list()
# numRuns<-1
# for(c in 0:(numConfs-1)){
#   print(paste0("run-",((c+21)*numRuns):(((c+21)*numRuns)+numRuns-1)))
#   runs <- c(runs,list(paste0("run-",((c+21)*numRuns):(((c+21)*numRuns)+numRuns-1))))
# }
# i <- 1
# allMs <- data.table()
# for(nms in numMS){
#   print(paste(nms,"MS"))
#   r <- runs[[i]]
#   celllist <- list()
#   for(cl in 0:(numCells-1)){
#     celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
#       return(fread(lrun)[Delay>0.41])
#     }))
#   }
#   allCells <- rbindlist(celllist)
#   allCells[,numMs:=as.factor(nms)]
#   allCells[,MS:=NULL]
#   ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1.0),stddev=sd(Delay),period="10ms"),by=c("numMs")]))
#   i <- i+1
#   allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
#   allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
#   allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
# }
# pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
#   geom_step(direction="hv")+
#   ylab("Pr[ packet delay < X ]")+
#   xlab("X in ms")+
#   theme_bw()
# jpeg(filename=paste0(figdir,figname,"_10ms.jpeg"),quality=100,type="cairo")
# plot(pl)
# dev.off()
# cairo_pdf(filename=paste0(figdir,figname,"_10ms.pdf"))
# plot(pl)
# dev.off()
# 
# write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)