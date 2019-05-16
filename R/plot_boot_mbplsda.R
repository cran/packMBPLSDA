
# graphical function for boot mbplsda

###########################################################################

plot_boot_mbplsda <- function(obj, filename="PlotBootstrapMbplsda", propbestvar=0.5){
  
  appel     <- as.list(obj$call)

  pdf(paste0(filename,".pdf"), paper="a4r", width=12, height=12)
  
  # points and margins
  formespch <- rep(21:25,3)
  bgpch     <- rep(c(0,0.5,1),5)
  bgpch     <- sapply(bgpch,gray, alpha=1)
  
  par(mai=c(1,1,1,1))
  
  # 1. BIPc : means and IC
  if(nrow(obj$bipc)>1){
    obj$bipc <- obj$bipc
    blocsBIP <- factor(obj$bipc$blocks)
    plot(obj$bipc$mean, xlab="blocks", ylab="BIPc", main="Means and 95% CI of cumulated BIP", xaxt="n", 
         pch=formespch[as.numeric(blocsBIP)], bg=bgpch[as.numeric(blocsBIP)],
         las=2, ylim=c(min(obj$bipc[,"95CIinf"],na.rm=TRUE),max(obj$bipc[,"95CIsup"],na.rm=TRUE)))
    axis(1, at=seq(from=1, to=nrow(obj$bipc)), labels=(obj$bipc$blocks),
         cex.axis=min(1,(50/nrow(obj$bipc))))
    for(i in 1: nrow(obj$bipc)) {segments(x0=i,y0=obj$bipc[i,"95CIinf"], x1=i ,y1=obj$bipc[i,"95CIsup"], lty=as.numeric(blocsBIP)[i])}
  }
   
  # 2. VIPc
  ## all
  hist(obj$vipc$mean, freq=F, main="Histogram and density curve of mean cumulated VIP", xlab="cumulated VIP")
  abline(v=mean(obj$vipc$mean, na.rm=T), lwd=4)
  lines(density(obj$vipc$mean), lwd=2)
  legend("topright",c("density curve","mean cumulated VIP"), lwd=c(2,4), cex=0.8)
  
  ## higher than the 1-propbestvar quantile 
  obj$vipc <- obj$vipc[order(obj$vipc$mean, decreasing=FALSE),]
  blocsVIP <- factor(obj$vipc$block)
  vipcBest <- obj$vipc[obj$vipc$mean>quantile(obj$vipc$mean, probs=(1-propbestvar), na.rm=TRUE),]
  blocsVIPBest <- blocsVIP[obj$vipc$mean>quantile(obj$vipc$mean, probs=(1-propbestvar), na.rm=TRUE)]
  
  ### histogram
  if(is.null(length(vipcBest$mean))==FALSE){
    hist(vipcBest$mean, freq=F, main=paste("Histogram and density curve \n of the",(propbestvar)*100,"% best mean cumulated VIP"), xlab="cumulated VIP")
    lines(density(vipcBest$mean), lwd=2)
  }
  
  ### means and CI
  margebas=min(11,max(6,(800/nrow(vipcBest))))
  par(mar=c(margebas, 4, 4, 2))
  plot(vipcBest$mean, xlab="", ylab="VIPc", main=paste("Means and 95% CI of the",(propbestvar)*100,"% best mean cumulated VIP"), 
       xaxt="n", pch=formespch[as.numeric(blocsVIPBest)], bg=bgpch[as.numeric(blocsVIPBest)], cex=0.8, las=2,
       ylim=c(min(vipcBest[,"95CIinf"],na.rm=TRUE),max(vipcBest[,"95CIsup"],na.rm=TRUE)))
  axis(1, at=seq(from=1, to=nrow(vipcBest)), labels=vipcBest$variables, las=2, 
       cex.axis=min(0.8,(50/nrow(vipcBest))) )
  mtext("variables", side=1, line=round((margebas-2),0))
  for(i in 1: nrow(vipcBest)) {segments(x0=i,y0=vipcBest[i,"95CIinf"], x1=i ,y1=vipcBest[i,"95CIsup"], lty=as.numeric(blocsVIPBest)[i])}
  legend("bottomright",paste("block",levels(blocsVIP)), cex=0.8, lty=as.numeric(factor(levels(blocsVIP))), 
         pch = formespch[as.numeric(factor(levels(blocsVIP)))],pt.bg = bgpch[as.numeric(factor(levels(blocsVIP)))] )
  par(mar=c(5, 4, 4, 2))
  
  # 3. loadings = faX
  
  for(i in 1:length(obj$faX)){
    matfaX       <- obj$faX[[i]]
    matfaX       <- matfaX[order(matfaX$mean, decreasing=FALSE),]
    blocsfaX     <- factor(matfaX$block)
    faXBest      <- matfaX[abs(matfaX$mean)>quantile(abs(matfaX$mean), probs=(1-propbestvar), na.rm=TRUE),]
    blocsfaXBest <- blocsfaX[abs(matfaX$mean)>quantile(abs(matfaX$mean), probs=(1-propbestvar), na.rm=TRUE)]
    
    ## faX histogram with all variables
    hist(matfaX$mean, freq=F, main=paste("Histogram and density curve of variables loadings on axe",i), xlab="loadings")
    lines(density(matfaX$mean), lwd=2)
    
    ## means and CI of variables with the higher loadings (absolute values)
    margebas=min(11,max(6,(800/nrow(faXBest))))
    par(mar=c(margebas, 4, 4, 2))
    plot(faXBest$mean, xlab="", ylab="loadings", main=paste("Mean and 95% CI of the",(propbestvar)*100,"% best variables loadings on axe",i), xaxt="n", 
         pch=formespch[as.numeric(blocsfaXBest)], bg=bgpch[as.numeric(blocsfaXBest)],cex=0.8, las=2,
         ylim=c(min(faXBest[,"95CIinf"],na.rm=TRUE),max(faXBest[,"95CIsup"],na.rm=TRUE)))
    axis(1, at=seq(from=1, to=nrow(faXBest)), labels=faXBest$variables, las=2, 
         cex.axis=min(0.8,(50/nrow(faXBest))) )
    mtext("variables", side=1, line=(margebas-2))
    for(i in 1: nrow(faXBest)) {segments(x0=i,y0=faXBest[i,"95CIinf"], x1=i ,y1=faXBest[i,"95CIsup"], lty=as.numeric(blocsfaXBest)[i])}
    legend("bottomright",paste("block",levels(blocsfaX)), cex=0.8, lty=as.numeric(factor(levels(blocsfaX))), 
           pch = formespch[as.numeric(factor(levels(blocsfaX)))], pt.bg = bgpch[as.numeric(factor(levels(blocsfaX)))] )
    
    par(mar=c(5, 4, 4, 2))
  }
  ## loadings plot
  matfaXall                         <- matrix(NA, nrow=nrow(obj$faX[[1]]), ncol=1+length(obj$faX), dimnames = list((obj$faX[[1]])$variables, c("blocks",names(obj$faX))))
  blocsfaXall                       <- factor(obj$faX[[1]]$block)
  matfaXall[,1:(length(obj$faX)+1)] <- cbind(blocsfaXall,sapply(1:length(obj$faX), function(x) obj$faX[[x]]$mean))
  
  ### IF MORE THAN ONE DIMENSION
  if(length(obj$faX)>1){
    if(length(obj$faX)%%2!=0) {matfaXall <- cbind(matfaXall[,1:length(obj$faX)],matfaXall[,length(obj$faX):(length(obj$faX)+1)])}
    
    minifaX <- min(matfaXall[,2:(dim(matfaXall)[2])])
    maxifaX <- max(matfaXall[,2:(dim(matfaXall)[2])])
    
    for (i in seq(from=2, to=dim(matfaXall)[2], by=2)){
      plot(matfaXall[,i],matfaXall[,i+1], main="Loadings plot (mean values)",
           pch=formespch[as.numeric(matfaXall[,"blocks"])],  bg=bgpch[as.numeric(matfaXall[,"blocks"])],
           xlab=colnames(matfaXall)[i], ylab=colnames(matfaXall)[i+1], cex=0.8, 
           xlim=c(minifaX,maxifaX)*1.05, ylim=c(minifaX,maxifaX)*1.05,
           las=1)
      abline(h=0,v=0)
      legend("bottomright",paste("block",levels(blocsfaXall)), cex=0.8, 
             pch = formespch[as.numeric(factor(levels(blocsfaXall)))], pt.bg=bgpch[as.numeric(factor(levels(blocsfaXall)))])
    }
  }
  ### IF ONE DIMENSION
  if(length(obj$faX)==1){
    minifaX <- min(matfaXall[,2])
    maxifaX <- max(matfaXall[,2])
    
    plot(matfaXall[,2],main="Loadings plot (mean values)",
         pch=formespch[as.numeric(matfaXall[,"blocks"])],  bg=bgpch[as.numeric(matfaXall[,"blocks"])], 
         xlab="variables", ylab=colnames(matfaXall)[2], cex=0.8, 
         ylim=c(minifaX,maxifaX)*1.05,
         las=1)
    abline(h=0)
    legend("bottomright",paste("block",levels(blocsfaXall)), cex=0.8, 
           pch = formespch[as.numeric(factor(levels(blocsfaXall)))], pt.bg = bgpch[as.numeric(factor(levels(blocsfaXall)))])
    
  }
  
  # 4. regression coefficients = XYcoef
  
  for(i in 1:length(obj$XYcoef)){
    matXYcoef       <- obj$XYcoef[[i]]
    matXYcoef       <- matXYcoef[order(matXYcoef$mean, decreasing=FALSE),]
    blocsXYcoef     <- factor(matXYcoef$block)
    XYcoefBest      <- matXYcoef[abs(matXYcoef$mean)>quantile(abs(matXYcoef$mean), probs=(1-propbestvar), na.rm=TRUE),]
    blocsXYcoefBest <- blocsXYcoef[abs(matXYcoef$mean)>quantile(abs(matXYcoef$mean), probs=(1-propbestvar), na.rm=TRUE)]
    
    ## XYcoef histogram with all variables
    hist(matXYcoef$mean, freq=F, main=paste("Histogram and density curve of variables regression coefficients \n to predict category",names(obj$XYcoef)[i]), xlab="regression coefficients")
    lines(density(matXYcoef$mean), lwd=2)
    
    ## means and CI of variables with the higher regression coefficients (absolute values)
    margebas=min(11,max(6,(800/nrow(XYcoefBest))))
    par(mar=c(margebas, 4, 4, 2))
    plot(XYcoefBest$mean, xlab="", ylab="regression coefficients", main=paste("Mean and 95% CI of the",(propbestvar)*100,"% best variables regression coefficients \n to predict category",names(obj$XYcoef)[i]), 
         xaxt="n", pch=formespch[as.numeric(blocsXYcoefBest)], bg=bgpch[as.numeric(blocsXYcoefBest)],cex=0.8, las=2,
         ylim=c(min(XYcoefBest[,"95CIinf"],na.rm=TRUE),max(XYcoefBest[,"95CIsup"],na.rm=TRUE)))
    axis(1, at=seq(from=1, to=nrow(XYcoefBest)), labels=XYcoefBest$variables, las=2, 
         cex.axis=min(0.8,(50/nrow(XYcoefBest))) )
    mtext("variables", side=1, line=(margebas-2))
    for(i in 1: nrow(XYcoefBest)) {segments(x0=i,y0=XYcoefBest[i,"95CIinf"], x1=i ,y1=XYcoefBest[i,"95CIsup"], lty=as.numeric(blocsXYcoefBest)[i])}
    legend("bottomright",paste("block",levels(blocsXYcoef)), cex=0.8, lty=as.numeric(factor(levels(blocsXYcoef))), 
           pch = formespch[as.numeric(factor(levels(blocsXYcoef)))], pt.bg = bgpch[as.numeric(factor(levels(blocsXYcoef)))])
    
    par(mar=c(5, 4, 4, 2))
  }
  
  dev.off()
}