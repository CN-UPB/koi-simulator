#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# It expects to find the KoI scripts/ directory in it's parent directory.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

source("../scripts/kbest_delays.R")

numCells <- 9

delays_ms_exploration(c(5,10,15,20),30,"./results_full/",9,"delays_5_to_20_ms.pdf","delays_smaller_1ms_5_to_20_ms.tab")
delays_ms_exploration(c(20,30,40,50),1,"./results_k_1/",9,"delays_k_1.pdf","delays_smaller_1ms_k_1.tab")

delays_k_exploration <- function(){
  numMS <- c(20,30,40,50)
  numRuns <- 1
  resDir <- "./results_k_exploration/"
  numCells <- 9
  kVals <- c(1,2,3,4,5)
  figname <- "delays_k_exploration.pdf"
  tablename <- "delays_smaller_1ms_by_k.tab"
  numConfs <- length(numMS)*length(kVals)
  figdir <- paste0(resDir,"../figures/")
  runs <- list()
  for(c in 0:(numConfs-1)){
    print(c)
    print(paste0("run-",(c*numRuns):((c*numRuns)+numRuns-1)))
    runs <- c(runs,list(paste0("run-",(c*numRuns):((c*numRuns)+numRuns-1))))
  }
  ecdfs <- data.table()
  cairo_pdf(filename=paste0(figdir,figname),onefile = TRUE)
  i <- 1
  allCells <- data.table()
  for(nms in numMS){
    plotData <- data.table()
    print(paste(nms,"MS"))
    for(k in kVals){
      r <- runs[[i]]
      celllist <- list()
      for(cl in 0:(numCells-1)){
        celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
          return(fread(lrun)[Delay>0.41][,k:=k])
        }))
      }
      allCells <- rbindlist(celllist)
      allCells[,MS:=NULL]
      plotData <- rbind(plotData,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=k])
      plotData <- rbind(plotData,data.table(ecx=-Inf,ecy=0.0,k=k))
      plotData <- rbind(plotData,data.table(ecx=Inf,ecy=1.0,k=k))
      ecdfs <- rbind(ecdfs,unique(allCells[,list(numMs=as.factor(nms),lt1msDelay=ecdf(Delay)(1),stddev=sd(Delay)),by=c("k")]))
      i <- i+1
    }
    plotData[,k:=as.factor(k)]
    pl <- ggplot(plotData,aes(x=ecx,y=ecy,colour=k))+geom_step(direction="hv")+ylab("Pr[ packet delay < X ]")+xlab("X in ms")+
      ggtitle(paste0("Transmission Delay ECDFs ",nms," Mobile Stations"))
    plot(pl)
  }
  dev.off()
  write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)
}

delays_k_exploration()