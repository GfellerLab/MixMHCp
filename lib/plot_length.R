
cm <- commandArgs(trailingOnly = TRUE)

min <- as.numeric(cm[2])
max <- as.numeric(cm[3])
ncl_min <- as.numeric(cm[4])
ncl_max <- as.numeric(cm[5])
trash <- as.numeric(cm[6])
lb <- c(0:5)/5

for(i in ncl_min:ncl_max){
  
  m <- read.table(paste(cm[1],"weights/weights_", i,".txt", sep=""), header=T)

  for(j in 1:i){
    sum <- sum(m[,i+j+2+trash])
    if(sum>0){
      tm <- m[,i+j+2+trash]/sum
    } else {
      tm <- rep(0,length(m[,i+j+2+trash]))   
    }
    png(paste(cm[1],"weights/plots/lg_",i,"_",j,".png", sep=""))
    plot(m[,1], tm, xlab="", ylim=c(0,1), ylab="", cex.axis=1.5, yaxt="n", xaxt="n", type="b")
    axis(2, at=lb,labels=lb, las=1, cex.axis=1.5, cex=1.5)
    axis(1, at=min:max,labels=min:max, las=1, cex.axis=1.5, cex=1.5)
    title(paste("Motif ", j, ", N=",sum, sep=""), xlab="Peptide length", ylab="", cex.lab=2, cex.main=3)
    dev.off()
  }
  if(trash==1){
    sum <- sum(m[,2*i+3+trash])
    if(sum>0){
      tm <- m[,2*i+3+trash]/sum
    } else {
      tm <- rep(0,length(m[,i+j+2+trash]))   
    }
    png(paste(cm[1],"weights/plots/lg_",i,"_Trash.png", sep=""))
    plot(m[,1], tm, xlab="", ylim=c(0,1), ylab="", cex.axis=1.5, yaxt="n", xaxt="n", type="b")
    axis(2, at=lb,labels=lb, las=1, cex.axis=1.5, cex=1.5)
    axis(1, at=min:max,labels=min:max, las=1, cex.axis=1.5, cex=1.5)
    title(paste("Trash, N=",sum, sep=""), xlab="Peptide length", ylab="", cex.lab=2, cex.main=3)
    dev.off()
  }
}
