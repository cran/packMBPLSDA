# plot pred MBlocs
#########################

plot_pred_mbplsda <- function (obj, filename="PlotPredMbplsda", propbestvar=0.5){
  
  pdf(paste0(filename,".pdf"), paper="a4r", width=12, height=12)
  
  # pch and margins
  formespch <- rep(21:25,3)
  bgpch     <- rep(c(0,0.5,1),5)
  bgpch     <- sapply(bgpch,gray, alpha=1)
  par(mai=c(1,1,1,1))
  
  # parameters and arguments
  appel            <- as.list(obj$call)
  initialModel     <- eval.parent(appel$object)
  bloY             <- eval.parent(appel$bloY)
  q                <- sum(bloY)
  nNoBin           <- sum(bloY!=2)  # nombre de variables non binaires
  optdim           <- eval.parent(appel$optdim)
  threshold        <- eval.parent(appel$threshold)
  EigVal           <- eval.parent(initialModel)$eig
  if(is.null(eval.parent(appel$algo))==TRUE) (algo <- c("max","gravity","threshold"))
  if(is.null(eval.parent(appel$algo))==FALSE) (algo <- eval.parent(appel$algo))
  
 
  # 1. eig
  barplot(EigVal[1:optdim], las=0, xlab="Axes", ylab="Eigen values", main="Eigen values", names.arg=1:optdim)

  # 2. scatter plot of observations
  ## IF MORE THAN ONE DIMENSION
  if(optdim>1){
    if(optdim>2 & optdim%%2!=0) {
      obj$lX  <- cbind(obj$lX[,1:(dim(obj$lX)[2]-1)],obj$lX[,(dim(obj$lX)[2]-1):(dim(obj$lX)[2])])
    }
    
    miniScoresInd <- min(obj$lX[,(q+2):dim(obj$lX)[2]])
    maxiScoresInd <- max(obj$lX[,(q+2):dim(obj$lX)[2]])
    
    # scatter plot with coloration according to truth
    par(mfrow=c(2,2))
    for(j in seq(from=(q+2), to=dim(obj$lX)[2], by=2)){
      for(i in 1:(sum(bloY))){
        plot(obj$lX[,j],obj$lX[,j+1], 
             pch=c(16,16,1)[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))], 
             col=c("grey","black","black")[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))],
             main=paste("Observed scatterplot \n obs colored by",colnames(eval.parent(appel$object)$tabY)[i]), 
             xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab=colnames(obj$lX)[j], ylab=colnames(obj$lX)[j+1],
             las=1)
        legend("bottomright", legend = c("0","1"), pch=c(16,16), col=c("grey","black"))
        abline(h=0, v=0)
      }
    }
   
    # scatter plot with coloration according to PredY.max
    if("max" %in% algo){
      for(j in seq(from=(q+2), to=dim(obj$lX)[2], by=2)){
        for(i in (sum(bloY)+2):dim(obj$PredY.max)[2]){
          plot(obj$lX[,j],obj$lX[,j+1], pch=c(16,16,1)[factor(obj$PredY.max[,i])], col=c("grey","black","black")[factor(obj$PredY.max[,i])],
               main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.max)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(obj$lX)[j], ylab=colnames(obj$lX)[j+1],
               las=1)
          abline(h=0, v=0)
          legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("method=maximal value", cex=0.75)
        }
      }
    }
    
    # scatter plot with coloration according to PredY.gravity
    if("gravity" %in% algo){
      for(j in seq(from=(q+2), to=dim(obj$lX)[2], by=2)){
        for(i in (sum(bloY)+2):dim(obj$PredY.gravity)[2]){
          plot(obj$lX[,j],obj$lX[,j+1], pch=c(16,16,1)[factor(obj$PredY.gravity[,i])], col=c("grey","black","black")[factor(obj$PredY.gravity[,i])],
               main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.gravity)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(obj$lX)[j], ylab=colnames(obj$lX)[j+1],
               las=1)
          abline(h=0, v=0)
          legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("method=center of gravity", cex=0.75)
        }
      }
    }
    
    
    # scatter plot with coloration according to PredY.threshold
    if("threshold" %in% algo){
      for(j in seq(from=(q+2), to=dim(obj$lX)[2], by=2)){
        for(i in (sum(bloY)+2):dim(obj$PredY.threshold)[2]){
          plot(obj$lX[,j],obj$lX[,j+1], pch=c(16,16,1)[factor(obj$PredY.threshold[,i])], col=c("grey","black","black")[factor(obj$PredY.threshold[,i])],
               main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.threshold)[i]), 
               xlim=c(miniScoresInd,maxiScoresInd)*1.05, ylim=c(miniScoresInd,maxiScoresInd)*1.05,
               xlab=colnames(obj$lX)[j], ylab=colnames(obj$lX)[j+1],
               las=1)
          abline(h=0, v=0)
          legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
          mtext("method=threshold", cex=0.75)
        }
      }
    }
  }
  
  ## IF ONE DIMENSION
  if(optdim==1){
    miniScoresInd <- min(obj$lX[,(q+2)])
    maxiScoresInd <- max(obj$lX[,(q+2)])
    
    # scatter plot with coloration according to truth
    par(mfrow=c(2,2))
    for(i in 1:(sum(bloY))){
      plot(obj$lX[,(q+2)],
           pch=c(16,16,1)[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))], 
           col=c("grey","black","black")[as.numeric(factor(eval.parent(appel$object)$tabY[,i]))],
           main=paste("Observed scatterplot \n obs colored by",colnames(eval.parent(appel$object)$tabY)[i]), 
           ylim=c(miniScoresInd,maxiScoresInd)*1.05,
           xlab="observations", ylab=colnames(obj$lX)[(q+2)],
           las=1)
      legend("bottomright", legend = c("0","1"), pch=c(16,16), col=c("grey","black"))
      abline(h=0)
    }
    
    # scatter plot with coloration according to PredY.max
    if("max" %in% algo){
      for(i in (sum(bloY)+2):dim(obj$PredY.max)[2]){
        plot(obj$lX[,(q+2)], pch=c(16,16,1)[factor(obj$PredY.max[,i])], col=c("grey","black","black")[factor(obj$PredY.max[,i])],
             main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.max)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab=colnames(obj$lX)[(q+2)],
             las=1)
        abline(h=0)
        legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("method=maximal value", cex=0.75)
      }
    }
    
    
    # scatter plot with coloration according to PredY.gravity
    if("gravity" %in% algo){
      for(i in (sum(bloY)+2):dim(obj$PredY.gravity)[2]){
        plot(obj$lX[,(q+2)], pch=c(16,16,1)[factor(obj$PredY.gravity[,i])], col=c("grey","black","black")[factor(obj$PredY.gravity[,i])],
             main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.gravity)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab=colnames(obj$lX)[(q+2)],
             las=1)
        abline(h=0)
        legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("method=center of gravity", cex=0.75)
      }
    }
    
    
    # scatter plot with coloration according to PredY.threshold
    if("threshold" %in% algo){
      for(i in (sum(bloY)+2):dim(obj$PredY.threshold)[2]){
        plot(obj$lX[,(q+2)], pch=c(16,16,1)[factor(obj$PredY.threshold[,i])], col=c("grey","black","black")[factor(obj$PredY.threshold[,i])],
             main=paste("Observed scatterplot \n obs colored by",colnames(obj$PredY.threshold)[i]), 
             ylim=c(miniScoresInd,maxiScoresInd)*1.05,
             xlab="observations", ylab=colnames(obj$lX)[(q+2)],
             las=1)
        abline(h=0)
        legend("bottomright", legend = c("0","1","NA"), pch=c(16,16,1), col=c("grey","black","black"))
        mtext("method=threshold", cex=0.75)
      }
    }
  }
  
    
    # 3. BIPc
  
  par(mfrow=c(1,1))

  plot(obj$BIPc[,(optdim+1)], xlab="blocks", ylab="BIPc", main="Cumulated BIP", xaxt="n", 
       pch=formespch[as.numeric(factor(obj$BIPc$blocks))], bg=bgpch[as.numeric(factor(obj$BIPc$blocks))],
       las=2, ylim=c(max(0,min(obj$BIPc[,(optdim+1)],na.rm=TRUE)-0.1),min(max(obj$BIPc[,(optdim+1)],na.rm=TRUE)+0.1)))
  axis(1, at=seq(from=1, to=nrow(obj$BIPc)), labels=obj$BIPc$blocks,
       cex.axis=min(1,(50/nrow(obj$BIPc))))
  
  # 4. VIPc
  ## all
  hist(obj$VIPc[,optdim+2], freq=F, main="Histogram and density curve of cumulated VIP", xlab="cumulated VIP")
  abline(v=mean(obj$VIPc[,optdim+2], na.rm=TRUE), lwd=4)
  lines(density(obj$VIPc[,optdim+2]), lwd=2)
  legend("topright",legend = c("density curve","mean cumulated VIP"), lwd=c(2,4), cex=0.8)
  
  ## higher than the 1-propbestvar quantile 
  obj$VIPc <- obj$VIPc[order(obj$VIPc[,optdim+2], decreasing=FALSE),]
  blocsVIP <- factor(obj$VIPc$block)
  VIPcBest <- obj$VIPc[obj$VIPc[,optdim+2] > quantile(obj$VIPc[,optdim+2], probs=(1-propbestvar), na.rm=TRUE),]
  blocsVIPcBest <- blocsVIP[obj$VIPc[,optdim+2] > quantile(obj$VIPc[,optdim+2], probs=(1-propbestvar), na.rm=TRUE)]
  
  ### histogram
  if(is.null(length(VIPcBest[,optdim+2]))==FALSE){
    hist(VIPcBest[,optdim+2], freq=F, main=paste("Histogram and density curve \n of the",(propbestvar)*100,"% best cumulated VIP"), xlab="cumulated VIP")
    lines(density(VIPcBest[,optdim+2]), lwd=2)
  }
  
  
  ### values
  margebas=min(11,max(6,(800/nrow(VIPcBest))))
  par(mar=c(margebas, 4, 4, 2))
  plot(VIPcBest[,optdim+2], xlab="", ylab="VIPc", main=paste("The",(propbestvar)*100,"% best cumulated VIP"), xaxt="n", 
       pch=formespch[as.numeric(blocsVIPcBest)], bg=bgpch[as.numeric(blocsVIPcBest)], cex=0.8, las=2,
       ylim=c(min(VIPcBest[,optdim+2],na.rm=TRUE)*0.95,max(VIPcBest[,optdim+2],na.rm=TRUE)*1.05))
  axis(1, at=seq(from=1, to=nrow(VIPcBest)), labels=VIPcBest$variables, las=2, 
       cex.axis=min(0.8,(50/nrow(VIPcBest))) )
  mtext("variables", side=1, line=(margebas-2))
  legend("bottomright",legend = paste("block",levels(blocsVIP)), cex=0.8, 
         pch = formespch[as.numeric(factor(levels(blocsVIP)))], pt.bg= bgpch[as.numeric(factor(levels(blocsVIP)))])
  par(mar=c(5, 4, 4, 2))
  
  
  # 4. loadings 
  
  for(i in 3:(optdim+2)){
    matfaX       <- obj$faX[,c(1,2,i)]
    matfaX       <- matfaX[order(matfaX[,3], decreasing=FALSE),]
    blocsfaX     <- factor(matfaX$block)
    faXBest      <- matfaX[abs(matfaX[,3]) > quantile(abs(matfaX[,3]), probs=(1-propbestvar), na.rm=TRUE),]
    blocsfaXBest <- blocsfaX[abs(matfaX[,3]) > quantile(abs(matfaX[,3]), probs=(1-propbestvar), na.rm=TRUE)]
    
    ## faX histogram with all variables
    hist(matfaX[,3], freq=F, main=paste("Histogram and density curve of variables loadings on axe",(i-2)), xlab="loadings")
    lines(density(matfaX[,3]), lwd=2)
    
    ## means and CI of half variables which have the highest loadings in absolute values
    margebas=min(11,max(6,(800/nrow(faXBest))))
    par(mar=c(margebas, 4, 4, 2))
    plot(faXBest[,3], xlab="", ylab="loadings", main=paste("The",(propbestvar)*100,"% best variables loadings on axe",(i-2)), xaxt="n", 
         pch=formespch[as.numeric(blocsfaXBest)], bg=bgpch[as.numeric(blocsfaXBest)], cex=0.8, las=2,
         ylim=c(min(faXBest[,3],na.rm=TRUE)*1.05,max(faXBest[,3],na.rm=TRUE)*1.05))
    axis(1, at=seq(from=1, to=nrow(faXBest)), labels=faXBest$variables, las=2, 
         cex.axis=min(0.8,(50/nrow(faXBest))) )
    mtext("variables", side=1, line=(margebas-2))
    legend("bottomright",legend = paste("block",levels(blocsfaX)), cex=0.8, 
           pch = formespch[as.numeric(factor(levels(blocsfaX)))], pt.bg = bgpch[as.numeric(factor(levels(blocsfaX)))])
    
    par(mar=c(5, 4, 4, 2))
  }
  
  ## loadings plots with pch according to blocks
  # IF MORE THAN ONE DIMENSION
  if (optdim >1){
    if(optdim%%2!=0 & optdim>2) {
      obj$faX <- cbind(obj$faX[,1:(dim(obj$faX)[2]-1)],obj$faX[,(dim(obj$faX)[2]-1):(dim(obj$faX)[2])])
    }
    
    minifaX       <- min(obj$faX[,3:dim(obj$faX)[2]])
    maxifaX       <- max(obj$faX[,3:dim(obj$faX)[2]])
    
    for(j in seq(from=3, to=dim(obj$faX)[2], by=2)){
      plot(obj$faX[,j],obj$faX[,j+1], main=paste("Loadings plot"), 
           pch=formespch[as.numeric(factor(obj$faX$block))], bg=bgpch[as.numeric(factor(obj$faX$block))], cex=0.8,
           xlim=c(minifaX,maxifaX)*1.05, ylim=c(minifaX,maxifaX)*1.05,
           xlab=colnames(obj$faX)[j], ylab=colnames(obj$faX)[j+1],
           las=1)
      legend("bottomright", legend = levels(obj$faX$block), cex=0.8, 
             pch=formespch[as.numeric(factor(levels(factor(obj$faX$block))))], pt.bg=bgpch[as.numeric(factor(levels(factor(obj$faX$block))))])
      abline(h=0, v=0)
    }
  }
  
  # IF ONE DIMENSION
  if (optdim == 1){
    minifaX       <- min(obj$faX[,3])
    maxifaX       <- max(obj$faX[,3])
    
    plot(obj$faX[,3],cex=0.8,main=paste("Loadings plot"),
         pch=formespch[as.numeric(factor(obj$faX$block))], bg=bgpch[as.numeric(factor(obj$faX$block))],
         ylim=c(minifaX,maxifaX)*1.05,
         xlab="variables", ylab=colnames(obj$faX)[3],
         las=1)
    legend("bottomright", legend = levels(obj$faX$block), cex=0.8, 
           pch=formespch[as.numeric(factor(levels(factor(obj$faX$block))))],pt.bg=bgpch[as.numeric(factor(levels(factor(obj$faX$block))))])
    abline(h=0)
  }
  
  
  
  # 5. regression coefficients = XYcoef
  
  for(i in 3:(q+2)){
    matXYcoef       <- obj$XYcoef[,c(1,2,i)]
    matXYcoef       <- matXYcoef[order(matXYcoef[,3], decreasing=FALSE),]
    blocsXYcoef     <- factor(matXYcoef$block)
    XYcoefBest      <- matXYcoef[abs(matXYcoef[,3]) > quantile(abs(matXYcoef[,3]), probs=(1-propbestvar), na.rm=TRUE),]
    blocsXYcoefBest <- blocsXYcoef[abs(matXYcoef[,3]) > quantile(abs(matXYcoef[,3]), probs=(1-propbestvar), na.rm=TRUE)]
    
    ## XYcoef histogram with all variables
    hist(matXYcoef[,3], freq=F, main=paste("Histogram and density curve of variables regression coefficients \n to predict category",
                                           colnames(obj$XYcoef)[i]), xlab="regression coefficients")
    lines(density(matXYcoef[,3]), lwd=2)
    
    ## means and IC of half variables which have the highest XYcoef in absolute values
    margebas=min(11,max(6,(800/nrow(XYcoefBest))))
    par(mar=c(margebas, 4, 4, 2))
    plot(XYcoefBest[,3], xlab="", ylab="regression coefficients", main=paste("The",(propbestvar)*100,"% best variables regression coefficients \n to predict category",names(obj$XYcoef)[i]), 
         xaxt="n", cex=0.8, las=2,
         pch=formespch[as.numeric(blocsXYcoefBest)], bg=bgpch[as.numeric(blocsXYcoefBest)], 
         ylim=c(min(XYcoefBest[,3],na.rm=TRUE)*1.05,max(XYcoefBest[,3],na.rm=TRUE)*1.05))
    axis(1, at=seq(from=1, to=nrow(XYcoefBest)), labels=XYcoefBest$variables, las=2, 
         cex.axis=min(0.8,(50/nrow(XYcoefBest))))
    mtext("variables", side=1, line=(margebas-2))
    legend("bottomright",legend = paste("block",levels(blocsXYcoef)), cex=0.8, 
           pch = formespch[as.numeric(factor(levels(blocsXYcoef)))], pt.bg=bgpch[as.numeric(factor(levels(blocsXYcoef)))])
    
    par(mar=c(5, 4, 4, 2))
  }
  
  dev.off()
}