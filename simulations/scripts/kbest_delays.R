require(data.table)
require(ggplot2)

reorganize_per_ms_files <- function(numMS,numRuns,resDir,numCells){
  numConfs <- length(numMS)
  runs <- list()
  for(c in 0:(numConfs-1)){
    runs <- c(runs,list(paste0("run-",(c*numRuns):((c*numRuns)+numRuns-1))))
  }
  i <- 1
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    for(cl in 0:(numCells-1)){
      for(run in r){
        msList <- data.table()
        for(msId in 0:(nms-1)){
          f <- paste0(resDir,"delay_","ms-",cl,"-",msId,"_",run,".dat")
          tmp <- fread(f)
          tmp[,MS:=msId]
          msList <- rbind(msList,tmp)
        }
        msList[,Delay:=V1][,V1:=NULL]
        write.table(msList,paste0(resDir,"delays-cell-",cl,"_",run,".dat"),sep="\t",row.names = FALSE)
      }
    }
    i <- i+1
  }
}

delays_ms_exploration <- function(numMS,numRuns,resDir,numCells,
                                  figname="delays_per_num_ms.pdf",
                                  tablename="delays_smaller_1ms_by_ms.tab"){
  numConfs <- length(numMS)
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
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1),stddev=sd(Delay)),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
    geom_step(direction="hv")+
    ylab("Pr[ packet delay < X ]")+
    xlab("X in ms")+
    ggtitle(paste0("Transmission Delay ECDFs"))
  plot(pl)
  dev.off()
  write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)
}

delays_period_exploration <- function(numMS,numRuns,resDir,numCells,
                                  figname="delays_per_num_ms.pdf",
                                  tablename="delays_smaller_1ms_by_ms.tab"){
  figdir <- paste0(resDir,"../figures/")
  periods = c("2ms","5ms")
  ecdfs <- data.table()
  numConfs <- length(numMS)
  # 2 ms Runs
  runs <- list()
  for(c in 0:(numConfs-1)){
    print(paste0("run-",(c*numRuns):((c*numRuns)+numRuns-1)))
    runs <- c(runs,list(paste0("run-",(c*numRuns):((c*numRuns)+numRuns-1))))
  }
  i <- 1
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1.0),stddev=sd(Delay),period="2ms"),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
    geom_step(direction="hv")+
    ylab("Pr[ packet delay < X ]")+
    xlab("X in ms")+
    theme_bw()
  jpeg(filename=paste0(figdir,figname,"_2ms.jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,"_2ms.pdf"))
  plot(pl)
  dev.off()

  # 5 ms Runs  
  runs <- list()
  for(c in 0:(numConfs-1)){
    print(paste0("run-",((c+7)*numRuns):(((c+7)*numRuns)+numRuns-1)))
    runs <- c(runs,list(paste0("run-",((c+7)*numRuns):(((c+7)*numRuns)+numRuns-1))))
  }
  i <- 1
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1.0),stddev=sd(Delay),period="5ms"),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
    geom_step(direction="hv")+
    ylab("Pr[ packet delay < X ]")+
    xlab("X in ms")+
    theme_bw()
  jpeg(filename=paste0(figdir,figname,"_5ms.jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,"_5ms.pdf"))
  plot(pl)
  dev.off()
  
  # 1 ms Runs
  runs <- list()
  for(c in 0:(numConfs-1)){
    print(paste0("run-",((c+14)*numRuns):(((c+14)*numRuns)+numRuns-1)))
    runs <- c(runs,list(paste0("run-",((c+14)*numRuns):(((c+14)*numRuns)+numRuns-1))))
  }
  i <- 1
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1.0),stddev=sd(Delay),period="1ms"),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
    geom_step(direction="hv")+
    ylab("Pr[ packet delay < X ]")+
    xlab("X in ms")+
    theme_bw()
  jpeg(filename=paste0(figdir,figname,"_1ms.jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,"_1ms.pdf"))
  plot(pl)
  dev.off()
  
  # 10 ms Runs
  runs <- list()
  for(c in 0:(numConfs-1)){
    print(paste0("run-",((c+21)*numRuns):(((c+21)*numRuns)+numRuns-1)))
    runs <- c(runs,list(paste0("run-",((c+21)*numRuns):(((c+21)*numRuns)+numRuns-1))))
  }
  i <- 1
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1.0),stddev=sd(Delay),period="10ms"),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+
    geom_step(direction="hv")+
    ylab("Pr[ packet delay < X ]")+
    xlab("X in ms")+
    theme_bw()
  jpeg(filename=paste0(figdir,figname,"_10ms.jpeg"),quality=100,type="cairo")
  plot(pl)
  dev.off()
  cairo_pdf(filename=paste0(figdir,figname,"_10ms.pdf"))
  plot(pl)
  dev.off()
  
  write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)
}

delays_ms_multifloor <- function(numMS,numRuns,resDir,numCells,cellsPerFloor,
                                  figname="delays_per_num_ms.pdf",
                                  tablename="delays_smaller_1ms_by_ms.tab"){
  numConfs <- length(numMS)
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
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1),stddev=sd(Delay)),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+geom_step(direction="hv")+ylab("Pr[ packet delay < X ]")+xlab("X in ms")+
    ggtitle(paste0("Transmission Delay ECDFs"))
  plot(pl)
  dev.off()
  write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)
}

deadlines_passed <- function(figname,tablename,title){
  figdir <- "./figures/"
  jpeg(filename=paste0(figdir,figname),quality=100,type="cairo")
  periods = c("1ms","2ms","5ms","10ms")
  data <- fread(paste0(figdir,tablename))
  pl <- ggplot(data,aes(x=numMs,y=lt1msDelay,colour=period))+
    geom_line(aes(group=period))+
    geom_point()+
    ylab("Pr[ packet delay < 1ms ]")+
    xlab("Nr. of Devices")+
    theme_linedraw()
  plot(pl)
  
  dev.off()
}

delays_ms_multifloor <- function(numMS,numRuns,resDir,numCells,cellsPerFloor,
                                 figname="delays_per_num_ms.pdf",
                                 tablename="delays_smaller_1ms_by_ms.tab"){
  numConfs <- length(numMS)
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
  allMs <- data.table()
  for(nms in numMS){
    print(paste(nms,"MS"))
    r <- runs[[i]]
    celllist <- list()
    for(cl in 0:(numCells-1)){
      celllist <- c(celllist,lapply(paste0(resDir,"delays-","cell-",cl,"_",runs[[i]],".dat"), function(lrun){
        return(fread(lrun)[Delay>0.41])
      }))
    }
    allCells <- rbindlist(celllist)
    allCells[,numMs:=as.factor(nms)]
    allCells[,MS:=NULL]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(lt1msDelay=ecdf(Delay)(1),stddev=sd(Delay)),by=c("numMs")]))
    i <- i+1
    allMs <- rbind(allMs,allCells[,list(ecx=unique(Delay),ecy=ecdf(Delay)(unique(Delay))),by=numMs])
    allMs <- rbind(allMs,data.table(ecx=-Inf,ecy=0.0,numMs=nms))
    allMs <- rbind(allMs,data.table(ecx=Inf,ecy=1.0,numMs=nms))
  }
  pl <- ggplot(allMs,aes(x=ecx,y=ecy,colour=numMs))+geom_step(direction="hv")+ylab("Pr[ packet delay < X ]")+xlab("X in ms")+
    ggtitle(paste0("Transmission Delay ECDFs"))
  plot(pl)
  dev.off()
  write.table(ecdfs,paste0(figdir,tablename),sep="\t",row.names = FALSE)
}