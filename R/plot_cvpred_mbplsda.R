
# plot CV pred MBlocs
#########################

plot_cvpred_mbplsda <- function (obj, filename="PlotCVpredMbplsda"){
  
  pdf(paste0(filename,".pdf"), paper="a4r", width=12, height=12)
  par(mai=c(1,1,1,1))
  
  # parameters and arguments
  appel            <- as.list(obj$call)
  initialModel     <- eval.parent(appel$object)
  bloY             <- eval.parent(appel$bloY)
  nNoBin           <- sum(bloY!=2)  # nombre de variables non binaires
  optdim           <- eval.parent(appel$optdim)
  if(is.null(eval.parent(appel$algo))==TRUE) (algo <- c("max","gravity","threshold"))
  if(is.null(eval.parent(appel$algo))==FALSE) (algo <- eval.parent(appel$algo))
  initialScoresInd <- eval.parent(initialModel)$lX[,1:optdim]
  
  
  # IF MORE THAN ONE DIMENSION
  if(optdim >1){
    if(optdim%%2!=0 & optdim>2) {initialScoresInd <- cbind(initialScoresInd[,1:(optdim-1)],initialScoresInd[,(optdim-1):(optdim)])}
    
    miniScoresInd <- min(initialScoresInd)
    maxiScoresInd <- max(initialScoresInd)
    
    # scatter plot with coloration according to the true statut 
    par(mfrow=c(2,2))
    for(j in seq(from=1, to=dim(initialScoresInd)[2], by=2)){
      for(i in 1:(sum(bloY))){
        plot(initialScoresInd[,j],initialScoresInd[,j+1], 
             pch=c(16,16,1)[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))], 
             col=c("grey","black","black")[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))],
             main=paste("Observed scatterplot \n obs colored by",colnames(eval.parent(appel$object)$tabY)[i]), 
             xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab=colnames(initialScoresInd)[j], ylab=colnames(initialScoresInd)[j+1],
             las=1)
        abline(h=0,v=0)
        legend("bottomright", c("0","1"), pch=c(16,16), col=c("grey","black"))
      }
    }
    
    # scatter plot with coloration according to matPredYv.max
    if("max" %in% algo){
      for(j in seq(from=1, to=dim(initialScoresInd)[2], by=2)){
        for(i in (sum(bloY)+2):(dim(obj$matPredYv.max)[2])){
          plot(initialScoresInd[,j],initialScoresInd[,j+1], pch=c(16,16,1)[factor(obj$matPredYv.max[,i])], col=c("grey","black","black")[factor(obj$matPredYv.max[,i])],
               main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.max)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(initialScoresInd)[j], ylab=colnames(initialScoresInd)[j+1],
               las=1)
          abline(h=0,v=0)
          legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("subset=validation, method=maximal value", cex=0.75)
        }
      }
    }
    
    
    # scatter plot with coloration according to matPredYv.gravity
    if("gravity" %in% algo){
      for(j in seq(from=1, to=dim(initialScoresInd)[2], by=2)){
        for(i in (sum(bloY)+2):(dim(obj$matPredYv.gravity)[2])){
          plot(initialScoresInd[,j],initialScoresInd[,j+1], pch=c(16,16,1)[factor(obj$matPredYv.gravity[,i])], col=c("grey","black","black")[factor(obj$matPredYv.gravity[,i])],
               main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.gravity)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(initialScoresInd)[j], ylab=colnames(initialScoresInd)[j+1],
               las=1)
          abline(h=0,v=0)
          legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("subset=validation, method=center of gravity", cex=0.75)
        }
      }
    }
    
    
    if("threshold" %in% algo){
      # scatter plot with coloration according to matPredYv.threshold
      for(j in seq(from=1, to=dim(initialScoresInd)[2], by=2)){
        for(i in (sum(bloY)+2):(dim(obj$matPredYv.threshold)[2])){
          plot(initialScoresInd[,j],initialScoresInd[,j+1], pch=c(16,16,1)[factor(obj$matPredYv.threshold[,i])], col=c("grey","black","black")[factor(obj$matPredYv.threshold[,i])],
               main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.threshold)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(initialScoresInd)[j], ylab=colnames(initialScoresInd)[j+1],
               las=1)
          abline(h=0,v=0)
          legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("subset=validation, method=threshold", cex=0.75)
        }
      }
    }
  }
  
  
  # IF ONE DIMENSION
  if(optdim==1){
    miniScoresInd <- min(initialScoresInd)
    maxiScoresInd <- max(initialScoresInd)
    
    # scatter plot with coloration according to the true statut
    par(mfrow=c(2,2))
    
    for(i in 1:(sum(bloY))){
      plot(initialScoresInd, 
           pch=c(16,16,1)[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))], 
           col=c("grey","black","black")[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))],
           main=paste("Observed scatterplot \n obs colored by",colnames(eval.parent(appel$object)$tabY)[i]), 
           ylim=c(miniScoresInd,maxiScoresInd)*1.05,
           xlab="observations", ylab="Ax1",
           las=1)
      abline(h=0)
      legend("bottomright", c("0","1"), pch=c(16,16), col=c("grey","black"))
    }
    
    # scatter plot with coloration according to matPredYv.max
    if("max" %in% algo){
      for(i in (sum(bloY)+2):(dim(obj$matPredYv.max)[2])){
        plot(initialScoresInd, pch=c(16,16,1)[factor(obj$matPredYv.max[,i])], col=c("grey","black","black")[factor(obj$matPredYv.max[,i])],
             main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.max)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab="Ax1",
             las=1)
        abline(h=0)
        legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("subset=validation, method=maximal value", cex=0.75)
      }
    }
    
    
    # scatter plot with coloration according to matPredYv.gravity
    if("gravity" %in% algo){
      for(i in (sum(bloY)+2):(dim(obj$matPredYv.gravity)[2])){
        plot(initialScoresInd, pch=c(16,16,1)[factor(obj$matPredYv.gravity[,i])], col=c("grey","black","black")[factor(obj$matPredYv.gravity[,i])],
             main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.gravity)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab="Ax1",
             las=1)
        abline(h=0)
        legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("subset=validation, method=center of gravity", cex=0.75)
      }
    }
    
    if("threshold" %in% algo){
      # scatter plot with coloration according to matPredYv.threshold
      for(i in (sum(bloY)+2):(dim(obj$matPredYv.threshold)[2])){
        plot(initialScoresInd, pch=c(16,16,1)[factor(obj$matPredYv.threshold[,i])], col=c("grey","black","black")[factor(obj$matPredYv.threshold[,i])],
             main=paste("Observed scatterplot obs colored \n by cross validated",colnames(obj$matPredYv.threshold)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab="Ax1",
             las=1)
        abline(h=0)
        legend("bottomright", c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("subset=validation, method=threshold", cex=0.75)
      }
    }
  }
  
  dev.off()
  
}