
cm <- commandArgs(trailingOnly = TRUE)

min <- as.numeric(cm[2])
max <- as.numeric(cm[3])
ncl <- as.numeric(cm[4])
trash <- as.numeric(cm[5])
lb <- c(0:5)/5

for(i in 1:ncl){
  
  m <- read.table(paste(cm[1],"weights/weights_", i,".txt", sep=""), header=T)

  for(j in 1:i){
    png(paste(cm[1],"weights/plots/lg_",i,"_",j,".png", sep=""))
    plot(m[,1], m[,i+j+2+trash], xlab="", ylim=c(0,1), ylab="", cex.axis=1.5, yaxt="n", xaxt="n", type="b")
    axis(2, at=lb,labels=lb, las=1, cex.axis=1.5, cex=1.5)
    axis(1, at=min:max,labels=min:max, las=1, cex.axis=1.5, cex=1.5)
    title(paste("Motif ", j, sep=""), xlab="Peptide length", ylab="", cex.lab=2, cex.main=3)
    dev.off()
  }
  

}