
# graphical fonction for MBPLSDA permutation test

###########################################################################

plot_permut_mbplsda <- function(obj, filename="PlotPermutationTest", MainPlot="Permutation test results \n (subset of validation)"){
  
  appel     <- as.list(obj$call)
  threshold <- eval.parent(appel$threshold)
  npermut   <- eval.parent(appel$npermut)
  dimlabP   <- c("NoPermut",paste("permut", (1 : npermut), sep = ""))
  bloY      <- eval.parent(appel$bloY)
  if(is.null(eval.parent(appel$algo))==TRUE) (algo <- c("max","gravity","threshold"))
  if(is.null(eval.parent(appel$algo))==FALSE) (algo <- eval.parent(appel$algo))
  nNoBin    <- sum(bloY!=2)  # nb no binary variables
  if(is.null(eval.parent(appel$outputs))==TRUE) (outputs <- c("ER","ConfMat","AUC"))
  if(is.null(eval.parent(appel$outputs))==FALSE) (outputs <- eval.parent(appel$outputs))
  if((("ER" %in% outputs) == FALSE) & (("AUC" %in% outputs) == FALSE)) 
    stop("'ER' or 'AUC' expected in outputs of permut_mbplsda")
  
  pdf(paste0(filename,".pdf"), width=12, height=12)
  par(mai=c(1,1,1,1))
  par(mfrow=c(2,2))
  
  if("ER" %in% outputs){
    # plot ERvMax
    if("max" %in% algo){
      regressionM           <- lm(obj$ErrorRateVglobal.max[(obj$ErrorRateVglobal.max$variable=="global" & obj$ErrorRateVglobal.max$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
                                  ~ obj$prctGlob.Ychange.values[2:(npermut+1),2])
      
      plot(obj$ErrorRateVglobal.max[(obj$ErrorRateVglobal.max$variable=="global" & obj$ErrorRateVglobal.max$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
           ~ obj$prctGlob.Ychange.values[2:(npermut+1),2], col="grey", xlim=c(0,1), ylim=c(0,1),pch=16,
           xlab="% modified values in the Y-block", ylab="mean error rate of classification", main=MainPlot, las=1, cex=0.8)
      abline(regressionM, lty=1, lwd=1)
      points(obj$prctGlob.Ychange.values[1,2],obj$ErrorRateVglobal.max[(obj$ErrorRateVglobal.max$variable=="global" & obj$ErrorRateVglobal.max$dataP == dimlabP[1]),"mean"], 
             col="black", pch=16)
      legend("bottomright", c("ERv.max","ERv.max regression line","ERv.max without permut."), cex=1, 
             lty=c(NA,1,NA), lwd=c(NA,1,NA), 
             pch=c(16,NA,16),col=c("grey","black","black"), pt.cex=c(0.8,NA,1))
    }
    
    # plot ERvGravity
    if("gravity" %in% algo){
      regressionG           <- lm(obj$ErrorRateVglobal.gravity[(obj$ErrorRateVglobal.gravity$variable=="global" & obj$ErrorRateVglobal.gravity$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
                                  ~ obj$prctGlob.Ychange.values[2:(npermut+1),2])
      
      plot(obj$ErrorRateVglobal.gravity[(obj$ErrorRateVglobal.gravity$variable=="global" & obj$ErrorRateVglobal.gravity$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
           ~ obj$prctGlob.Ychange.values[2:(npermut+1),2], col="grey", xlim=c(0,1), ylim=c(0,1),pch=16,las=1,
           xlab="% modified values in the Y-block Y", ylab="mean error rate of classification", main=MainPlot, las=1, cex=0.8)
      abline(regressionG, lty=1, lwd=1)
      points(obj$prctGlob.Ychange.values[1,2],obj$ErrorRateVglobal.gravity[(obj$ErrorRateVglobal.gravity$variable=="global" & obj$ErrorRateVglobal.gravity$dataP == dimlabP[1]),"mean"], 
             col="black", pch=16)
      legend("bottomright", c("ERv.gravity","ERv.gravity regression line","ERv.gravity without permut."), cex=1, 
             lty=c(NA,1,NA), lwd=c(NA,1,NA), 
             pch=c(16,NA,16),col=c("grey","black","black"), pt.cex=c(0.8,NA,1))
    }
     
    # plot ERvThreshold
    if("threshold" %in% algo){
      regressionT           <- lm(obj$ErrorRateVglobal.threshold[(obj$ErrorRateVglobal.threshold$variable=="global" & obj$ErrorRateVglobal.threshold$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
                                  ~ obj$prctGlob.Ychange.values[2:(npermut+1),2])
      
      plot(obj$ErrorRateVglobal.threshold[(obj$ErrorRateVglobal.threshold$variable=="global" & obj$ErrorRateVglobal.threshold$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
           ~ obj$prctGlob.Ychange.values[2:(npermut+1),2], col="grey", xlim=c(0,1), ylim=c(0,1),pch=16,las=1,
           xlab="% modified values in the Y-block Y", ylab="mean error rate of classification", main=MainPlot, las=1, cex=0.8)
      abline(regressionT, lty=1, lwd=1)
      points(obj$prctGlob.Ychange.values[1,2],obj$ErrorRateVglobal.threshold[(obj$ErrorRateVglobal.threshold$variable=="global" & obj$ErrorRateVglobal.threshold$dataP == dimlabP[1]),"mean"], 
             col="black", pch=16)
      legend("bottomright", c("ERv.threshold","ERv.threshold regression line","ERv.threshold without permut."),
             pch=c(16,NA,16), pt.cex=c(0.8,NA,1), 
             lty=c(NA,1,NA), lwd=c(NA,1,NA), cex=1, 
             col=c("grey","black","black"))
    }
  }
  
  # plot AUC
  if((nNoBin==0) & ("AUC" %in% outputs)){
    regressionAUC         <- lm(obj$AUCv.global[(obj$AUCv.global$variable=="Mean" & obj$AUCv.global$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
                                ~ obj$prctGlob.Ychange.values[2:(npermut+1),2])

    plot(obj$AUCv.global[(obj$AUCv.global$variable=="Mean" & obj$AUCv.global$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]
         ~ obj$prctGlob.Ychange.values[2:(npermut+1),2], col="grey", xlim=c(0,1), ylim=c(0,1),pch=16,las=1,
         xlab="% modified values in the Y-block Y", ylab="mean AUC", main=MainPlot, las=1, cex=0.8)
    abline(regressionAUC, lty=1, lwd=1)
    points(obj$prctGlob.Ychange.values[1,2],obj$AUCv.global[(obj$AUCv.global$variable=="Mean" & obj$AUCv.global$dataP == dimlabP[1]),"mean"], 
           col="black", pch=16, cex=1)
    legend("topright", c("AUCv","AUCv regression line","AUCv without permut."),
           pch=c(16,NA,16), pt.cex=c(0.8,NA,1), 
           lty=c(NA,1,NA), lwd=c(NA,1,NA), cex=1, 
           col=c("grey","black","black"))
  }
   
  dev.off()
}