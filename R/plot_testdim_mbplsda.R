
# graphic of results of testdim mbplsda

plot_testdim_mbplsda <- function (obj, filename="PlotTestdimMbplsda"){
  
  appel     <- as.list(obj$call)
  threshold <- eval.parent(appel$threshold)
  bloY      <- eval.parent(appel$bloY)
  nNoBin    <- sum(bloY!=2)  # nb no binary variables
  if(is.null(eval.parent(appel$algo))==TRUE) (algo <- c("max","gravity","threshold"))
  if(is.null(eval.parent(appel$algo))==FALSE) (algo <- eval.parent(appel$algo))
  nf        <- eval.parent(appel$object)$nf
  
  if(is.null(eval.parent(appel$outputs))==TRUE) (outputs <- c("ER","ConfMat","AUC"))
  if(is.null(eval.parent(appel$outputs))==FALSE) (outputs <- eval.parent(appel$outputs))
  if((("ER" %in% outputs) == FALSE) & (("AUC" %in% outputs) == FALSE)) 
    stop("'ER' or 'AUC' expected in outputs of testdim_mbplsda")

  # means
  if("ER" %in% outputs){
    ERvM.mean <- data.frame(globalERvM=obj$ErrorRateVglobal.max[obj$ErrorRateVglobal.max$variable=="global","mean"]) 
    ERcM.mean <- data.frame(globalERcM=obj$ErrorRateCglobal.max[obj$ErrorRateCglobal.max$variable=="global","mean"])
    ERvG.mean <- data.frame(globalERvG=obj$ErrorRateVglobal.gravity[obj$ErrorRateVglobal.gravity$variable=="global","mean"])
    ERcG.mean <- data.frame(globalERcG=obj$ErrorRateCglobal.gravity[obj$ErrorRateCglobal.gravity$variable=="global","mean"])
    ERvT.mean <- data.frame(globalERvT=obj$ErrorRateVglobal.threshold[obj$ErrorRateVglobal.threshold$variable=="global","mean"]) 
    ERcT.mean <- data.frame(globalERcT=obj$ErrorRateCglobal.threshold[obj$ErrorRateCglobal.threshold$variable=="global","mean"])
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    AUCv.mean <- data.frame(meanAUCv=obj$AUCv.global[obj$AUCv.global$variable=="Mean","mean"])
    AUCc.mean <- data.frame(meanAUCc=obj$AUCc.global[obj$AUCc.global$variable=="Mean","mean"])
  }

  erreurs <- data.frame(dimlab=paste0(rep("Ax",nf),1:nf))
  
  if("ER" %in% outputs){
    if("max" %in% algo){
      erreurs <- cbind(erreurs,ERvM.mean, ERcM.mean)
    }
    if("gravity" %in% algo){
      erreurs <- cbind(erreurs,ERvG.mean, ERcG.mean)
    }
    if("threshold" %in% algo){
      erreurs <- cbind(erreurs,ERvT.mean, ERcT.mean)
    }
  }
  
  if((nNoBin==0) & ("AUC" %in% outputs)){
    auc <- cbind(AUCv.mean, AUCc.mean)
  }
  
  # CI inf
  if("ER" %in% outputs){
    ERvM.ICinf <- data.frame(globalERvM=obj$ErrorRateVglobal.max[obj$ErrorRateVglobal.max$variable=="global","95CIinf"]) 
    ERcM.ICinf <- data.frame(globalERcM=obj$ErrorRateCglobal.max[obj$ErrorRateCglobal.max$variable=="global","95CIinf"])
    ERvG.ICinf <- data.frame(globalERvG=obj$ErrorRateVglobal.gravity[obj$ErrorRateVglobal.gravity$variable=="global","95CIinf"])
    ERcG.ICinf <- data.frame(globalERcG=obj$ErrorRateCglobal.gravity[obj$ErrorRateCglobal.gravity$variable=="global","95CIinf"])
    ERvT.ICinf <- data.frame(globalERvT=obj$ErrorRateVglobal.threshold[obj$ErrorRateVglobal.threshold$variable=="global","95CIinf"]) 
    ERcT.ICinf <- data.frame(globalERcT=obj$ErrorRateCglobal.threshold[obj$ErrorRateCglobal.threshold$variable=="global","95CIinf"])
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    AUCv.ICinf <- data.frame(meanAUCv=obj$AUCv.global[obj$AUCv.global$variable=="Mean","95CIinf"])
    AUCc.ICinf <- data.frame(meanAUCc=obj$AUCc.global[obj$AUCc.global$variable=="Mean","95CIinf"])
  }
  
  erreursICinf <- data.frame(dimlab=paste0(rep("Ax",nf),1:nf))
  if("ER" %in% outputs){
    if("max" %in% algo){
      erreursICinf <- cbind(erreursICinf,ERvM.ICinf, ERcM.ICinf)
    }
    if("gravity" %in% algo){
      erreursICinf <- cbind(erreursICinf,ERvG.ICinf, ERcG.ICinf)
    }
    if("threshold" %in% algo){
      erreursICinf <- cbind(erreursICinf,ERvT.ICinf, ERcT.ICinf)
    }
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    aucICinf <- cbind(AUCv.ICinf, AUCc.ICinf)
  }
 
  # CI sup
  if("ER" %in% outputs){
    ERvM.ICsup <- data.frame(globalERvM=obj$ErrorRateVglobal.max[obj$ErrorRateVglobal.max$variable=="global","95CIsup"]) 
    ERcM.ICsup <- data.frame(globalERcM=obj$ErrorRateCglobal.max[obj$ErrorRateCglobal.max$variable=="global","95CIsup"])
    ERvG.ICsup <- data.frame(globalERvG=obj$ErrorRateVglobal.gravity[obj$ErrorRateVglobal.gravity$variable=="global","95CIsup"])
    ERcG.ICsup <- data.frame(globalERcG=obj$ErrorRateCglobal.gravity[obj$ErrorRateCglobal.gravity$variable=="global","95CIsup"])
    ERvT.ICsup <- data.frame(globalERvT=obj$ErrorRateVglobal.threshold[obj$ErrorRateVglobal.threshold$variable=="global","95CIsup"]) 
    ERcT.ICsup <- data.frame(globalERcT=obj$ErrorRateCglobal.threshold[obj$ErrorRateCglobal.threshold$variable=="global","95CIsup"])
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    AUCv.ICsup <- data.frame(meanAUCv=obj$AUCv.global[obj$AUCv.global$variable=="Mean","95CIsup"])
    AUCc.ICsup <- data.frame(meanAUCc=obj$AUCc.global[obj$AUCc.global$variable=="Mean","95CIsup"])
  }
  
  erreursICsup <- data.frame(dimlab=paste0(rep("Ax",nf),1:nf))
  if("ER" %in% outputs){
    if("max" %in% algo){
      erreursICsup <- cbind(erreursICsup,ERvM.ICsup, ERcM.ICsup)
    }
    if("gravity" %in% algo){
      erreursICsup <- cbind(erreursICsup,ERvG.ICsup, ERcG.ICsup)
    }
    if("threshold" %in% algo){
      erreursICsup <- cbind(erreursICsup,ERvT.ICsup, ERcT.ICsup)
    }
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    aucICsup <- cbind(AUCv.ICsup, AUCc.ICsup)
  }
  
  # breaks
  if("ER" %in% outputs){
    ERmin  <- min(erreursICinf[,-1], na.rm=TRUE)*0.9
    ERmax  <- max(erreursICsup[,-1], na.rm=TRUE)*1.1
  }
  
  if((nNoBin==0) & ("AUC" %in% outputs)){
    AUCmin <- min(aucICinf, na.rm=TRUE)*0.9
    AUCmax <- max(aucICsup, na.rm=TRUE)*1.1
  }
  
  
  # graph
  
  pdf(paste0(filename,".pdf"), paper = "a4r", width=12, height=12)
  
  par(mai=c(1,1,1,1))
  par(mfrow=c(2,2))
  
  if("ER" %in% outputs){
    
    if("max" %in% algo){
      plot(erreurs[,"globalERcM"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=calibration - method=maximal value", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERcM"], x1=i ,y1=erreursICsup[i,"globalERcM"], lty=1)}
      
      plot(erreurs[,"globalERvM"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=validation - method=maximal value", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERvM"], x1=i ,y1=erreursICsup[i,"globalERvM"], lty=1)}
    }
    
    if("gravity" %in% algo){
      plot(erreurs[,"globalERcG"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=calibration - method=center of gravity", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERcG"], x1=i ,y1=erreursICsup[i,"globalERcG"], lty=1)}
      
      plot(erreurs[,"globalERvG"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=validation - method=center of gravity", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERvG"], x1=i ,y1=erreursICsup[i,"globalERvG"], lty=1)}
    }
    
    if("threshold" %in% algo){
      plot(erreurs[,"globalERcT"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=calibration - method=threshold", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERcT"], x1=i ,y1=erreursICsup[i,"globalERcT"], lty=1)}
      
      plot(erreurs[,"globalERvT"], type="b", pch=16, lty=1, cex=0.8, ylim=c(ERmin,ERmax), xlab="number of components", ylab="error rate",xaxs="r", las=1,
           main="Mean and 95% CI of the global error rate \n subset=validation - method=threshold", cex.main=0.9, xaxp=c(1, nf, nf-1))
      for(i in 1: nf) {segments(x0=i,y0= erreursICinf[i,"globalERvT"], x1=i ,y1=erreursICsup[i,"globalERvT"], lty=1)}
    }
  }
  
  if((nNoBin==0) & ("AUC" %in% outputs)){
    plot(auc[,"meanAUCc"], type="b", pch=16, lty=1, cex=0.8,ylim=c(AUCmin,AUCmax), xlab="number of components", ylab="AUC",xaxs="r", las=1,
         main="Mean and 95% CI of the area under curve \n subset=calibration", cex.main=0.9, xaxp=c(1, nf, nf-1)) 
    for(i in 1: dim(auc)[1]) {segments(x0=i,y0= aucICinf[i,"meanAUCc"], x1=i ,y1=aucICsup[i,"meanAUCc"], lty=1)}
    
    plot(auc[,"meanAUCv"], type="b", pch=16, lty=1, cex=0.8,ylim=c(AUCmin,AUCmax), xlab="number of components", ylab="AUC",xaxs="r", las=1,
         main="Mean and 95% CI of the area under curve \n subset=validation", cex.main=0.9, xaxp=c(1, nf, nf-1)) 
    for(i in 1: dim(auc)[1]) {segments(x0=i,y0= aucICinf[i,"meanAUCv"], x1=i ,y1=aucICsup[i,"meanAUCv"], lty=1)}
  }

  dev.off()
}