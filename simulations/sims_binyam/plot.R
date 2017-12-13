#!/usr/bin/env Rscript
# This script will execute evaluation scripts for the local experiment.
# The script library can be found in https://git.cs.upb.de/koi/simulations.git.

require(data.table)
require(ggplot2)
require(plot3D)

assignRegion <- function(regions,sinr){
  return(regions[starts<sinr & ends>sinr]$region)
}

resDir <- "./results/"
figdir <- paste0(resDir,"../figures/")
figname <- paste0(figdir,"SINR.pdf")
numLevels <- 15

data <- fread(paste0(resDir,"coeff_table_down-0_run-0.dat"))[MS==0][RB==0][,TTI:=TTI*0.001]
bands <- seq(min(data$Coeff),max(data$Coeff),length.out=numLevels+1)
print(length(bands))
h <- hist(data$Coeff,breaks=bands)
h <- data.table(region=seq(0,length(bands)-2),starts=bands[0:15],ends=bands[2:16],count=h$counts)
data <- data[,list(region=assignRegion(h,Coeff),SINR=Coeff),by=TTI]
data[,nextRegion:=region[2:nrow(data)]]
print(data)
p <- ggplot(data,aes(x=TTI,y=SINR))+
  geom_line()+
  ylab("SNR [dB]")+
  xlab("Time [s]")+
  ggtitle("SNR over Time")+
  geom_hline(yintercept = seq(min(data$SINR),max(data$SINR),length.out=numLevels+1),colour="orange")+
  geom_hline(yintercept = min(data$SINR),colour="blue")+
  geom_hline(yintercept = max(data$SINR),colour="blue")
for(i in seq(1,nrow(h))){
  p <- p+annotate("text",x=0.0,y=(h$starts[i]+(h$ends[i]-h$starts[i])/2),label=paste0("State ",i-1),colour="red")
}
ph <- ggplot()+geom_bar(h,mapping = aes(x=region,y=count),stat="identity")+
  ylab("Number of Occurences")+xlab("State")+ggtitle("State Counts")
cairo_pdf(figname,onefile = TRUE)
plot(p)
plot(ph)
transitions <- data[,list(Probability=nrow(.SD)/nrow(data)),by=c("region","nextRegion")]
pl <- ggplot(transitions,aes(x=region,y=nextRegion,colour=Probability))+geom_point()+
  ylab("Current State")+xlab("Previous State")+ggtitle("State Transition Probabilities")+
  scale_x_continuous(breaks=c(0,seq(-1:numLevels+1)))+
  scale_y_continuous(breaks=c(0,seq(-1:numLevels+1)))
plot(pl)
dev.off()
jpeg(paste0(figdir,"snr_over_time.jpeg"))
plot(p)
dev.off()
jpeg(paste0(figdir,"state_counts.jpeg"))
plot(ph)
dev.off()
jpeg(paste0(figdir,"state_transitions.jpeg"))
plot(pl)
dev.off()
write.table(h,paste0(figdir,"counts.tab"),row.names = FALSE)