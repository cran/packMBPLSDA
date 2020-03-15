
# -----------------------------------------------------------------------------------------
# predictions for a given nb of components
# -----------------------------------------------------------------------------------------


pred_mbplsda <- function(object, optdim , threshold = 0.5, bloY, algo=c("max","gravity","threshold")){

# bloY = NEEDED VECTOR = nb of categories by Y block variable

#  library(pROC)
  
  if (!inherits(object, "mbplsda")) 
    stop("Object of type 'mbplsda' expected")
  
  # Arguments 
  appel  <- as.list(object$call)
  method <- as.character(appel[[1]])
  scale  <- eval.parent(appel$scale)
  option <- eval.parent(appel$option)
  if(class(try(eval.parent(appel$ktabX), silent = TRUE))=="try-error") {
    stop("ktabX must be in the Global Environment")
  }
  X      <- eval.parent(appel$ktabX)
  if(class(try(eval.parent(appel$dudiY), silent = TRUE))[1]=="try-error") {
    stop("dudiY must be in the Global Environment")
  }
  Y      <- eval.parent(appel$dudiY)  
  nr     <- nrow(Y$tab)  
  q      <- ncol(Y$tab)
  h      <- object$rank
  nblo   <- length(X$blo) 
  Ky     <- length(bloY)
  
  Var <- as.factor(rep(1 : Ky, bloY))
  #Mod <- unlist(sapply(1:Ky, function(x) rep(c(1:bloY[x]))))
  cnames <- colnames(Y$tab) #paste0("Var",Var,"Mod",Mod)
  
  blo    <- sapply(1:nblo, function(k) dim(X[[k]])[2]) # nb variables by X block
  Bloc   <- as.factor(rep(1:nblo, blo))
  
  rnamesX <- row.names(X)
  rnamesY <- row.names(Y$tab)
  
  # 1. outputs preparation 
  classY <- matrix(NA, nrow=nr, ncol=3*Ky)
  colnames(classY) <- c(paste0("Ymax_Var",1:Ky),paste0("Ycentroid_Var",1:Ky), paste0("Ythreshold_Var",1:Ky))   
  rownames(classY) <- rnamesY
    
  if("max" %in% algo){
    ClassPredY.max <- matrix(0,nr,ncol=q)
    colnames(ClassPredY.max) <- cnames
    rownames(ClassPredY.max) <- rnamesY
    
    AccuracyPredY.max <- matrix(0,nr,ncol=(q+Ky+1))
    colnames(AccuracyPredY.max) <- c(cnames,paste0("Var",1 : Ky),"global")
    rownames(AccuracyPredY.max) <- rnamesY
  }
  if("gravity" %in% algo){
    ClassPredY.gravity <- matrix(0,nr,ncol=q)
    colnames(ClassPredY.gravity) <- cnames
    rownames(ClassPredY.gravity) <- rnamesY
    
    AccuracyPredY.gravity <- matrix(0,nr,ncol=(q+Ky+1))
    colnames(AccuracyPredY.gravity) <- c(cnames,paste0("Var",1 : Ky),"global")
    rownames(AccuracyPredY.gravity) <- rnamesY
  }
  if("threshold" %in% algo){
    ClassPredY.threshold <- matrix(NA,nr,ncol=q)
    colnames(ClassPredY.threshold) <- cnames
    rownames(ClassPredY.threshold) <- rnamesY
    
    AccuracyPredY.threshold <- matrix(NA,nr,ncol=(q+Ky+1))
    colnames(AccuracyPredY.threshold) <- c(cnames,paste0("Var",1 : Ky),"global")
    rownames(AccuracyPredY.threshold) <- rnamesY
  }
  if(sum(bloY!=2)==0){
    AUC <- matrix(NA, ncol=q, nrow=3, dimnames=list(c("auc","infCI","supCI"), cnames))
    AUC.global <- matrix(NA, ncol=(Ky+1), nrow=3, dimnames=list(c("auc","infCI","supCI"), c(paste0("Var",1:Ky),"Mean")))
  }
  
  ConfMat.ErrorRate <- matrix(NA, nrow = 15, ncol= q, 
                              dimnames = list(c("TruePos.max","TrueNeg.max","FalsePos.max","FalseNeg.max","ErrorRate.max",
                                                "TruePos.gravity","TrueNeg.gravity","FalsePos.gravity","FalseNeg.gravity","ErrorRate.gravity",
                                                "TruePos.threshold","TrueNeg.threshold","FalsePos.threshold","FalseNeg.threshold","ErrorRate.threshold"), 
                                              cnames))

  #residY  <- matrix(NA, nrow = nr, ncol = q)
  #colnames(residY) <- cnames
  #rownames(residY) <- rnamesY
  
  ErrorRate.global <- matrix(NA, nrow=3, ncol=(Ky+1),
                             dimnames = list(c("ErrorRate.max","ErrorRate.gravity","ErrorRate.threshold"),
                                             c(paste0("Var",1:Ky),"global")))
  
  
  
  # 2. model
  # application of mbplsda    
  rescal <- mbplsda(dudiY = Y, ktabX = X, scale = scale, option = option, scannf = FALSE, nf = optdim)

  # coefficients for raw data
  XYcoef.cal    <- sapply(rescal$XYcoef.raw, function(x) x[, optdim])
  intercept.cal <- sapply(rescal$intercept.raw, function(x) x[, optdim])
  
  # raw matrix
  X.mat <- cbind.data.frame(lapply(unclass(X)[1:nblo], scale, center = FALSE, scale = FALSE), stringsAsFactors = TRUE)
  
   
  # nb observations by Y category
  nbY1 <- sapply(1:q, function(g) sum(Y$tab[, g] == 1))
  
  
  # 3. predictions Y and RMSE
  # predictions Y by optdim
  predY           <- matrix(rep(intercept.cal, each=nr), ncol=q) + as.matrix(X.mat) %*% XYcoef.cal
  colnames(predY) <- cnames

  
  # RMSE
  #residY       <- as.matrix(Y$tab) - predY 
  #RMSE         <- sqrt(sum(residY^2) / (nr * q))
  
  
  # 4. predictions and error rates by Y category with max
  if("max" %in% algo){
    for (k in 1:Ky){
      ## index 
      classY[,k] <- sapply(1:nr, function(n) which.max(predY[n,Var == k]))
      
      ## disjunctive table
      for(n in 1:nr){
        if(k==1) {ClassPredY.max[n,classY[n,1]] <- 1}
        if(k>1)  {ClassPredY.max[n,(sum(bloY[1:(k-1)])+classY[n,k])] <- 1}
      }
    }
    
    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.max[n , 1:q] <- sapply(1:q, function(Q)(1-((ClassPredY.max[n,Q]-Y$tab[n,Q])^2)))
    }
        
    ## confusion matrix and error rates
    ConfMat.ErrorRate ["TruePos.max",]     <- sapply(1:q, function(l) length(which(ClassPredY.max[,l]==1 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["TrueNeg.max",]     <- sapply(1:q, function(l) length(which(ClassPredY.max[,l]==0 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalsePos.max",]    <- sapply(1:q, function(l) length(which(ClassPredY.max[,l]==1 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalseNeg.max",]    <- sapply(1:q, function(l) length(which(ClassPredY.max[,l]==0 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["ErrorRate.max",]   <- sapply(1:q, function(l) (length(which(ClassPredY.max[,l]==1 & Y$tab[,l]==0)) + length(which(ClassPredY.max[,l]==0 & Y$tab[,l]==1)))/nr)
  }
  

  # 5. predictions and errors rates with gravity
  if("gravity" %in% algo){
    ## group barycenters on global components
    Gravity <- matrix(0, nrow=q, ncol=optdim, dimnames=list(cnames, 1:optdim))
    
    if (optdim == 1){
      Gravity[,1] <- (sapply(1:q, function(g){
        if(nbY1[g]>1) (mean(rescal$lX[Y$tab[, g] == 1, 1]))
        else if(nbY1[g]==1) (rescal$lX[Y$tab[, g] == 1, 1])
        else if(nbY1[g]==0) (NA)
      }))
    }else{
      Gravity[,1:optdim] <- (t(sapply(1:q, function(g) {
        if(nbY1[g]>1) (apply(rescal$lX[Y$tab[, g] == 1,1:optdim], 2, mean))
        else if(nbY1[g]==1) (rescal$lX[Y$tab[, g] == 1,1:optdim])
        else if(nbY1[g]==0) (rep(NA,times=optdim))
      })))
    }
    
    
    ## distances to barycenters
    dist.eucl.gravity.Y          <- matrix(NA, nrow=nr, ncol=q, dimnames=list(1:nr,1:q))
    dist.eucl.gravity.Y[,1:q]    <- sapply(1:q, function(g) apply((rescal$lX[,1:optdim] -  matrix(rep(Gravity[g, ], each = nr), nrow = nr))^2, 1, sum))
     
    
    for (k in 1:Ky){
      
      ## index of predicted categories
      classY[,(k+Ky)] <- sapply(1:nr, function(n) which.min(dist.eucl.gravity.Y[n,Var == k]))
            
      ## disjonctive table
      for(n in 1:nr){
        if(k==1) {ClassPredY.gravity[n,classY[n,(1+Ky)]] <- 1}
        if(k>1)  {ClassPredY.gravity[n,(sum(bloY[1:(k-1)])+classY[n,(k+Ky)])] <- 1}
      }
    }
    
    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.gravity[n , 1:q] <- sapply(1:q, function(Q)(1-((ClassPredY.gravity[n,Q]-Y$tab[n,Q])^2)))
    }
        
    ## confusion matrix and error rates
    ConfMat.ErrorRate ["TruePos.gravity",]     <- sapply(1:q, function(l) length(which(ClassPredY.gravity[,l]==1 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["TrueNeg.gravity",]     <- sapply(1:q, function(l) length(which(ClassPredY.gravity[,l]==0 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalsePos.gravity",]    <- sapply(1:q, function(l) length(which(ClassPredY.gravity[,l]==1 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalseNeg.gravity",]    <- sapply(1:q, function(l) length(which(ClassPredY.gravity[,l]==0 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["ErrorRate.gravity",]   <- sapply(1:q, function(l) (length(which(ClassPredY.gravity[,l]==1 & Y$tab[,l]==0)) + length(which(ClassPredY.gravity[,l]==0 & Y$tab[,l]==1)))/nr)
  }
  

  # 6. predictions and error rates with threshold
  if("threshold" %in% algo){

    ## predicted categories indicators
    ClassPredY.threshold[predY>=threshold] <- 1
    ClassPredY.threshold[predY<threshold]  <- 0

    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.threshold[n , 1:q] <- sapply(1:q, function(Q)(1-((ClassPredY.threshold[n,Q]-Y$tab[n,Q])^2)))
    }
        
    ## confusion matrix and error rates
    ConfMat.ErrorRate ["TruePos.threshold",]    <- sapply(1:q, function(l) length(which(ClassPredY.threshold[,l]==1 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["TrueNeg.threshold",]    <- sapply(1:q, function(l) length(which(ClassPredY.threshold[,l]==0 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalsePos.threshold",]   <- sapply(1:q, function(l) length(which(ClassPredY.threshold[,l]==1 & Y$tab[,l]==0))/nr)
    ConfMat.ErrorRate ["FalseNeg.threshold",]   <- sapply(1:q, function(l) length(which(ClassPredY.threshold[,l]==0 & Y$tab[,l]==1))/nr)
    ConfMat.ErrorRate ["ErrorRate.threshold",]  <- sapply(1:q, function(l) (length(which(ClassPredY.threshold[,l]==1 & Y$tab[,l]==0)) + length(which(ClassPredY.threshold[,l]==0 & Y$tab[,l]==1)))/nr)

    for (k in 1:Ky){
      ## index of predicted categories
      classY[,(2*Ky+k)] <- sapply(1:nr, function(n) if(sum(ClassPredY.threshold[n,Var == k], na.rm=T)==1) (which.max(ClassPredY.threshold[n,Var == k])))
    }
  }
    
    ## auc if all groups are in Y
  if(sum(bloY!=2)==0){  
    if(sum(nbY1==0)==0){
      AUC["auc",]   <- sapply(1:q, function(l) ci.auc(Y$tab[,l],predY[,l])[2])
      AUC["infCI",] <- sapply(1:q, function(l) ci.auc(Y$tab[,l],predY[,l])[1])
      AUC["supCI",] <- sapply(1:q, function(l) ci.auc(Y$tab[,l],predY[,l])[3])
    }
  }

  # 7. error rates by variable and overall
  
  if("max" %in% algo){
    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.max[n,(q+1):(q+Ky)] <- sapply(1:Ky, function(k)(min(1-(ClassPredY.max[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    AccuracyPredY.max[,"global"]       <- apply(AccuracyPredY.max[,1:(q+Ky)], min, MARGIN = 1)
    ## error rates by variable and overall
    ErrorRate.global["ErrorRate.max",] <- 1-(apply(AccuracyPredY.max[,(q+1):(q+Ky+1)], mean, MARGIN = 2))
  }
  if("gravity" %in% algo){
    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.gravity[n,(q+1):(q+Ky)] <- sapply(1:Ky, function(k)(min(1-(ClassPredY.gravity[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    AccuracyPredY.gravity[,"global"]       <- apply(AccuracyPredY.gravity[,1:(q+Ky)], min, MARGIN = 1)
    ## error rates by variable and overall
    ErrorRate.global["ErrorRate.gravity",] <- 1-(apply(AccuracyPredY.gravity[,(q+1):(q+Ky+1)], mean, MARGIN = 2))
  }
  if("threshold" %in% algo){
    ## accuracy indicators
    for(n in 1:nr){
      AccuracyPredY.threshold[n,(q+1):(q+Ky)] <- sapply(1:Ky, function(k)(min(1-(ClassPredY.threshold[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    AccuracyPredY.threshold[,"global"]       <- apply(AccuracyPredY.threshold[,1:(q+Ky)], min, MARGIN = 1)
    ## error rates by variable and overall
    ErrorRate.global["ErrorRate.threshold",] <- 1-(apply(AccuracyPredY.threshold[,(q+1):(q+Ky+1)], mean, MARGIN = 2))
  }
 
  ## AUC by variable and overall
  if(sum(bloY!=2)==0){
    AUC.global[c("auc","infCI","supCI"),(1:Ky)]<- AUC[c("auc","infCI","supCI"),seq(from=1, to=q, by=2)]
    if(Ky>1){
      testAUC                                  <- t.test(AUC.global["auc",(1:Ky)])
      AUC.global["auc","Mean"]                 <- testAUC$estimate 
      AUC.global["infCI","Mean"]               <- testAUC$conf.int[1]
      AUC.global["supCI","Mean"]               <- testAUC$conf.int[2]
    }else{
      AUC.global[,"Mean"]                      <- AUC.global[,1]
    }
  }  
  
    
  # 8. other results and results list
  predictions <- list()
  
  #predictions$RMSE = RMSE
  
  # coefficients, VIPc, BIPc, loadings, components
  
  block     <- unlist(list(sapply(1:nblo, function(b) rep(names(X$blo)[b], (X$blo)[b]))))
  variables <- unlist(list(sapply(1:nblo, function(v) colnames(X[[v]]))))
  
  predictions$XYcoef <- data.frame(variables,block,sapply(rescal$XYcoef, function(x) x[, optdim]), stringsAsFactors = TRUE)
  predictions$VIPc   <- data.frame(variables,block,rescal$vipc, stringsAsFactors = TRUE)
  predictions$BIPc   <- data.frame(blocks=names(X$blo),rescal$bipc, stringsAsFactors = TRUE)
  predictions$faX    <- data.frame(variables,block,rescal$faX, stringsAsFactors = TRUE)
  predictions$lX     <- data.frame(obs=rnamesY,Y$tab,rescal$lX, stringsAsFactors = TRUE)
  rownames(predictions$XYcoef) <- rownames(predictions$VIPc) <- rownames(predictions$BIPc) <- rownames(predictions$faX) <- rownames(predictions$lX) <- NULL
  
  # 9. results matrix
  ConfMat.ErrorRate                       <- ConfMat.ErrorRate[complete.cases(ConfMat.ErrorRate), ]
  predictions$ConfMat.ErrorRate           <- data.frame(rates=rownames(ConfMat.ErrorRate),ConfMat.ErrorRate, stringsAsFactors = TRUE)
  rownames(predictions$ConfMat.ErrorRate) <- NULL
  
  #ErrorRate.global                       <- ErrorRate.global[complete.cases(ErrorRate.global), ]
  predictions$ErrorRate.global           <- data.frame(rates=rownames(ErrorRate.global),ErrorRate.global, stringsAsFactors = TRUE)
  predictions$ErrorRate.global           <- predictions$ErrorRate.global[complete.cases(ErrorRate.global), ]
  rownames(predictions$ErrorRate.global) <- NULL
  
  if("max" %in% algo){
    predictions$PredY.max           <- cbind(rnamesY,Y$tab,ClassPredY.max,AccuracyPredY.max) 
    colnames(predictions$PredY.max) <- c("obs",paste0("Truth_",colnames(Y$tab)),paste0("Pred_",colnames(ClassPredY.max)),paste0("Accu_",colnames(AccuracyPredY.max)))
    rownames(predictions$PredY.max) <- NULL
  }
  if("gravity" %in% algo){
    predictions$PredY.gravity           <- cbind(rnamesY,Y$tab,ClassPredY.gravity,AccuracyPredY.gravity) 
    colnames(predictions$PredY.gravity) <- c("obs",paste0("Truth_",colnames(Y$tab)),paste0("Pred_",colnames(ClassPredY.gravity)),paste0("Accu_",colnames(AccuracyPredY.gravity)))
    rownames(predictions$PredY.gravity) <- NULL
  }
  if("threshold" %in% algo){
    predictions$PredY.threshold           <- cbind(rnamesY,Y$tab,ClassPredY.threshold,AccuracyPredY.threshold) 
    colnames(predictions$PredY.threshold) <- c("obs",paste0("Truth_",colnames(Y$tab)),paste0("Pred_",colnames(ClassPredY.threshold)),paste0("Accu_",colnames(AccuracyPredY.threshold)))
    rownames(predictions$PredY.threshold) <- NULL
  }

  if(sum(bloY!=2)==0){
    predictions$AUC           <- data.frame(values=rownames(AUC),AUC,AUC.global, stringsAsFactors = TRUE)
    rownames(predictions$AUC) <- NULL
  }

  predictions$call <- match.call()
  class(predictions) <- c("pred_mbplsda")
  predictions
}

