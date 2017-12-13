require(data.table)
require(ggplot2)

ratesOld <- function(){
  cairo_pdf(filename="Rates.pdf",onefile = TRUE)
  rdir <- "/home/michael/fgrn/koi/simulator/simulations/kbest_multi_line_scarcity/results/"
  runs <- c("run_0","run_1","run_2","run_3","run_4","run_5","run_6","run_7","run_8","run_9")
  #runs <- c("run_0")
  vec <- data.table()
  for(ms in c("ms_0_0","ms_0_1","ms_0_2","ms_0_3","ms_0_4","ms_0_5","ms_0_6","ms_0_7","ms_0_8","ms_0_9")){
    for(r in runs){
      fname <- paste0(rdir,r,"_rate_",ms,".dat")
      tmp <- fread(fname)
      tmp[,ms:=ms]
      vec <- rbind(vec,tmp)
    }
  }
  pl <- ggplot(vec,aes(V1,colour=ms))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X")+ggtitle("Transmission Rate ECDFs")
  plot(pl)
  dev.off()  
}

ratesNew <- function(){
  cairo_pdf(filename="Rates.pdf",onefile = TRUE)
  rdir <- "/home/michael/fgrn/koi/simulator/simulations/line_interference/results/"
  #runs <- c("run_0","run_1","run_2","run_3","run_4","run_5","run_6","run_7","run_8","run_9")
  runs <- c("run-0")
  vec <- data.table()
  for(ms in c("ms-0-0","ms-0-1","ms-0-2","ms-0-3","ms-0-4","ms-0-5","ms-0-6","ms-0-7","ms-0-8","ms-0-9")){
    for(r in runs){
      fname <- paste0(rdir,"rate-",ms,"_",r,".dat")
      tmp <- fread(fname)
      tmp[,ms:=ms]
      vec <- rbind(vec,tmp)
    }
  }
  pl <- ggplot(vec,aes(V1,colour=ms))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X")+ggtitle("Transmission Rate ECDFs")
  plot(pl)
  dev.off()  
}

rates_iteration <- function(){
  numConfs <- 3
  numMS <- 5
  rdir <- "/home/michael/fgrn/koi/simulator/simulations/shop_floor/results/"
  figdir <- paste0(rdir,"../figures/")
  runs <- list(list("run-0","run-1","run-2"),list("run-3","run-4","run-5"),list("run-6","run-7","run-8"))
  cells <- c(1,4,9)
  for(i in 1:numConfs){
    allCells <- data.table()
    ncells <- cells[i]
    cairo_pdf(filename=paste0(figdir,ncells,"_rates.pdf"),onefile = TRUE)
    for(cl in 0:(ncells-1)){
      vec <- data.table()
      for(ms in 0:(numMS-1)){
        for(r in runs[[i]]){
          fname <- paste0(rdir,"rate-","ms-",cl,"-",ms,"_",r,".dat")
          print(paste0("Conf:",i," Cell:",cl," MS:",ms," runs:",runs[[i]]))
          tmp <- fread(fname)
          tmp[,ms:=ms]
          vec <- rbind(vec,tmp)
        }
      }
      vec[,ms:=as.factor(ms)]
      pl <- ggplot(vec,aes(V1,colour=ms))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X")+
        ggtitle(paste0("Transmission Rate ECDFs ","Cell ",cl))
      plot(pl)
      vec[,ms:=NULL]
      vec[,cell:=as.factor(cl)]
      allCells <- rbind(allCells,vec)
    }
    pl <- ggplot(allCells,aes(V1,colour=cell))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X")+
      ggtitle(paste0("Transmission Rate ECDFs All Cells"))
    plot(pl)
    dev.off()
  }
}

rates_9x9 <- function(){
  numConfs <- 4
  numMS <- c(5,10,15,20)
  rdir <- "/home/michael/fgrn/koi/simulator/simulations/grid_9x9/results_full/"
  figdir <- paste0(rdir,"../figures/")
  runs <- list(paste0("run-",0:29),paste0("run-",30:59),paste0("run-",60:89),paste0("run-",90:119))
  cells <- 9
  ecdfs <- data.table()
  for(i in 1:numConfs){
    allCells <- data.table()
    nms <- numMS[i]
    cairo_pdf(filename=paste0(figdir,nms,"_rates.pdf"),onefile = TRUE)
    for(cl in 0:(cells-1)){
      vec <- data.table()
      for(ms in 0:(nms-1)){
        for(r in runs[[i]]){
          fname <- paste0(rdir,"rate-","ms-",cl,"-",ms,"_",r,".dat")
          print(paste0("Conf:",i," Cell:",cl," MS:",ms," runs:",runs[[i]]))
          tmp <- fread(fname)
          tmp[,inMS:=sum(V1),by= (seq(nrow(tmp)) - 1) %/% 4]
          vec <- rbind(vec,tmp)
        }
      }
      vec[,cell:=as.factor(cl)]
      allCells <- rbind(allCells,vec)
    }
    pl <- ggplot(allCells,aes(V1,colour=cell))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X per TTI")+
      ggtitle(paste0("Transmission Rate ECDFs All Cells"))
    plot(pl)
    pl <- ggplot(allCells,aes(inMS,colour=cell))+stat_ecdf(pad=FALSE)+ylab("Pr[ rate < x ]")+xlab("X per 1 ms")+
      ggtitle(paste0("Transmission Rate ECDFs All Cells"))
    plot(pl)
    dev.off()
    allCells[,numMs:=nms]
    ecdfs <- rbind(ecdfs,unique(allCells[,list(numMs=numMs,gt100bPerTTI=1-(ecdf(V1)(100)),gt100bPerMS=1-(ecdf(inMS)(100))),by=c("cell")]))
  }
  write.table(ecdfs,paste0(figdir,"rate_grater_100b.tab"),sep="\t",row.names = FALSE)
}
