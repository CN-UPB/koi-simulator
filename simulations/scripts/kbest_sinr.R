require(data.table)
require(ggplot2)
rdir <- "/home/michael/fgrn/koi/simulator/simulations/kbest_multi_line_scarcity/results/"
runs <- c("run_0","run_1","run_2","run_3","run_4","run_5","run_6","run_7","run_8","run_9")
#runs <- c("run_0")


sinrAll <- function(rdir,runs){
  cairo_pdf(filename="SINR_all.pdf",onefile = TRUE)
  vec <- data.table()
  for(ms in c("ms_0_0","ms_0_1","ms_0_2","ms_0_3","ms_0_4","ms_0_5","ms_0_6","ms_0_7","ms_0_8","ms_0_9")){
    for(r in runs){
      fname <- paste0(rdir,r,"_sinr_",ms,".dat")
      tmp <- fread(fname)
      tmp[,ms:=ms]
      vec <- rbind(vec,tmp[,list(SINR=SINR,ms=ms)])
    }
  }
  pl <- ggplot(vec,aes(SINR,colour=ms))+stat_ecdf(pad=FALSE)+ylab("Pr[ SINR < x ]")+xlab("X")+ggtitle("SINR ECDFs Over All RBs")
  plot(pl)
  dev.off()  
}

sinrScheduled <- function(rdir,runs){
  cairo_pdf(filename="SINR_scheduled.pdf",onefile = TRUE)
  vec <- data.table()
  for(cell in c("cell_0")){
    for(r in runs){
      fname <- paste0(rdir,r,"_schedule_up_",cell,".dat")
      tmp <- fread(fname)
      vec <- rbind(vec,tmp)
    }
  }
  vec[,MS:=factor(MS)]
  pl <- ggplot()+stat_ecdf(data=vec,aes(SINR,col=MS),pad=FALSE)+ylab("Pr[ SINR < x ]")+xlab("X")+ggtitle("SINR ECDFs Over Scheduled RBs")
  plot(pl)
  dev.off()  
}

sinrTogether <- function(rdir,runs){
  cairo_pdf(filename="SINR_together.pdf",onefile = TRUE)
  vec <- data.table()
  for(ms in c("ms_0_0","ms_0_1","ms_0_2","ms_0_3","ms_0_4","ms_0_5","ms_0_6","ms_0_7","ms_0_8","ms_0_9")){
    for(r in runs){
      fname <- paste0(rdir,r,"_sinr_",ms,".dat")
      tmp <- fread(fname)
      tmp[,ms:=ms]
      vec <- rbind(vec,tmp[,list(SINR=SINR,ms=ms)])
    }
  }
  vec[,t:="all"]
  for(cell in c("cell_0")){
    for(r in runs){
      fname <- paste0(rdir,r,"_schedule_up_",cell,".dat")
      tmp <- fread(fname)
      tmp <- tmp[,ms:=paste0("ms_0_",MS)][,t:="scheduled"][,list(ms,SINR,t)]
      vec <- rbind(vec,tmp)
    }
  }
  pl <- ggplot(vec,aes(SINR,colour=ms,linetype=t))+stat_ecdf(pad=FALSE)+ylab("Pr[ SINR < x ]")+xlab("X")+ggtitle("SINR ECDFs")
  plot(pl)
  dev.off()
}

sinr_all_iteration <- function(){
  numConfs <- 3
  numMS <- 5
  rdir <- "/home/michael/fgrn/koi/simulator/simulations/shop_floor/results/"
  figdir <- paste0(rdir,"../figures/")
  runs <- list(list("run-0","run-1","run-2"),list("run-3","run-4","run-5"),list("run-6","run-7","run-8"))
  cells <- c(1,4,9)
  for(i in 1:numConfs){
    allCells <- data.table()
    ncells <- cells[i]
    cairo_pdf(filename=paste0(figdir,ncells,"_sinr.pdf"),onefile = TRUE)
    for(cl in 0:(ncells-1)){
      vec <- data.table()
#      for(ms in 0:(numMS-1)){
#        for(r in runs[[i]]){
#          fname <- paste0(rdir,"sinr-","ms-",cl,"-",ms,"_",r,".dat")
#          print(paste0("Conf:",i," Cell:",cl," MS:",ms," runs:",runs[[i]]))
#          tmp <- fread(fname)
#          tmp[,ms:=ms]
#          vec <- rbind(vec,tmp[,list(SINR=SINR,ms=ms)])
#        }
#      }
#      vec[,t:="all"]
      for(r in runs[[i]]){
        print(paste0("Scheduled run ",r))
        fname <- paste0(rdir,"schedule-up-cell-",cl,"_",r,".dat")
        tmp <- fread(fname)
        tmp <- tmp[,ms:=MS][,t:="scheduled"][,list(ms,SINR,t)]
        vec <- rbind(vec,tmp)
      }
      vec[,ms:=as.factor(ms)]
      pl <- ggplot(vec,aes(SINR,colour=ms,linetype=t))+stat_ecdf(pad=FALSE)+ylab("Pr[ SINR < x ]")+xlab("X")+
        ggtitle(paste0("Scheduled SINR ECDFs"," Cell ",cl))
      plot(pl)
      vec[,ms:=NULL]
      vec[,cell:=as.factor(cl)]
      print("Binding to allCells")
      allCells <- rbind(allCells,vec)
    }
    pl <- ggplot(allCells,aes(SINR,colour=cell,linetype=t))+stat_ecdf(pad=FALSE)+ylab("Pr[ SINR < x ]")+xlab("X")+
      ggtitle(paste0("SINR ECDFs All Cells"," Cell ",cl))
    plot(pl)
    dev.off()
    #return(allCells)
  }
}

metis_test_boxplots <- function(){
  stations <- c("ms-0-0","ms-4-0")
  dat <- data.table()
  for(ms in stations){
    tmp <- fread(paste0("~/fgrn/koi/simulator/simulations/metis_performance/results/","sinr-",ms,"_run-0.dat"))[TTI<10]
    dat <- rbind(dat,tmp)
  }
  dat[,Cell:=as.factor(Cell)]
  dat[,TTI:=as.factor(TTI)]
  print(dat)
  pl <- ggplot()+geom_boxplot(data=dat,aes(TTI,SINR,color=Cell))
  plot(pl)
}

metis_test_rbs <- function(){
  stations <- c("ms-0-0")
  numRuns<-5
  runs<-paste0("run-",0:(numRuns-1))
  dat <- data.table()
  for(r in runs){
    for(ms in stations){
      tmp <- fread(paste0("~/fgrn/koi/simulator/simulations/metis_test_single_cell/results/","sinr-",ms,"_",r,".dat"))[TTI==3][,run:=r]
      dat <- rbind(dat,tmp)
    }    
  }
  pl <- ggplot()+geom_line(data=dat,aes(RB,SINR,color=run))+scale_color_discrete(labels=c("2 GHz","3 GHz","4 GHz","5 GHz","6 GHz"))
  plot(pl)
}

metis_test_rb_over_time <- function(){
  stations <- c("ms-0-0")
  dat <- data.table()
  for(ms in stations){
    tmp <- fread(paste0("~/fgrn/koi/simulator/simulations/metis_test_single_cell/results/","sinr-",ms,"_run-0.dat"))
    dat <- rbind(dat,tmp)
  }
  pl <- ggplot()+geom_line(data=dat[RB==50],aes(TTI,SINR))
  plot(pl)
}

metis_test_rbs()