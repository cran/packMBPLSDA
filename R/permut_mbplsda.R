
# -----------------------------------------------------------------------------------------
# permutation test
# -----------------------------------------------------------------------------------------


permut_mbplsda <- function(object, optdim, bloY, algo=c("max","gravity","threshold"), threshold = 0.5, nrepet = 100, npermut=100, nbObsPermut=NULL, outputs = c("ER","ConfMat","AUC"), cpus=1){

# nrepet = nb CV repetitions
# npermut = nb permutated Y blocks  
# nbObsPermut = nb of observations to permute to obtain each permutated Y block. If NULL, this nb is randomly chosen for each permutated Y block
# bloY = MUST BE GIVEN BY THE USER = vector of numbers of categories per variable in the Y block
  
  if (!inherits(object, "mbplsda")) 
    stop("Object of type'mbplsda' expected")
  
  # step 0. packages necessaires
  
  #library(microbenchmark) # to test time performances
  #library(parallel)  # for the repartition of jobs
  #library(doParallel) # for iterations
  #library(foreach)
  #library(pROC) # for auc
  #library(FactoMineR) # for RV coefficents
 
  
  # step 1. Arguments 
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
  Nc     <- round(2 * nr / 3)
  Nv     <- nr - Nc 
  
  Ky     <- length(bloY) # nb of Y variables
  
  Var    <- as.factor(rep(1 : Ky, bloY))
#  Mod    <- unlist(sapply(1:Ky, function(x) rep(c(1:bloY[x]))))
  cnames <- colnames(Y$tab) #paste0("Var",Var,"Mod",Mod)
  
  nblo   <- length(X$blo)  # nb X blocks
  blo    <- sapply(1:nblo, function(k) dim(X[[k]])[2]) # nb variables per X block
  Bloc   <- as.factor(rep(1:nblo, blo))
  
  nNoBin <- sum(bloY!=2)  # nb no binary variables 
  
  dimlabP <- c("NoPermut",paste("permut", (1 : npermut), sep = ""))
  
  # step 3. preparation of permutation outputs
  res                          <- list()
  res$RV.YYpermut.values       <- data.frame(dimlabP,rep(NA, (npermut+1)), stringsAsFactors = TRUE)
  colnames(res$RV.YYpermut.values) <- c("dimlabP","RV.YYpermut.values")
  
  res$cor.YYpermut.values      <- data.frame(dimlabP,matrix(NA, ncol=q, nrow=(npermut+1)), stringsAsFactors = TRUE) #, dimnames=list(dimlabP, cnames))
  colnames(res$cor.YYpermut.values) <- c("dimlabP",cnames)
  
  res$prctGlob.Ychange.values <- data.frame(dimlabP,rep(NA, (npermut+1)), stringsAsFactors = TRUE)
  colnames(res$prctGlob.Ychange.values) <- c("dimlabP","prctGlob.Ychange.values")
  
  res$prct.Ychange.values     <- data.frame(dimlabP, matrix(NA, ncol=q, nrow=(npermut+1)), stringsAsFactors = TRUE) #, dimnames=list(dimlabP, cnames))
  colnames(res$prct.Ychange.values) <- c("dimlabP",cnames)
  
  ## means
  if("RMSE" %in% outputs){
  	RMSECglobal.P        <- list() 
  	RMSEC.P              <- list() 
  	RMSEVglobal.P        <- list() 
  	RMSEV.P              <- list() 
  }
  
  if("max" %in% algo){
    if("ConfMat" %in% outputs){
      TruePosC.max <- TruePosV.max <- TrueNegC.max <- TrueNegV.max <- FalsePosC.max <- FalsePosV.max <- FalseNegC.max <- FalseNegV.max <- list()
    }
    if("ER" %in% outputs){
      ErrorRateC.max <- ErrorRateV.max <- list()
      ErrorRateCglobal.max <- ErrorRateVglobal.max <- list()
    }
  }
  
  if("gravity" %in% algo){
    if("ConfMat" %in% outputs){
      TruePosC.gravity <- TruePosV.gravity <- TrueNegC.gravity <- TrueNegV.gravity <- FalsePosC.gravity <- FalsePosV.gravity <- FalseNegC.gravity <- FalseNegV.gravity <- list()
    }
    if("ER" %in% outputs){
      ErrorRateC.gravity <- ErrorRateV.gravity <- list()
      ErrorRateCglobal.gravity <- ErrorRateVglobal.gravity <- list()
    }
  }
  
  if("threshold" %in% algo){
    if("ConfMat" %in% outputs){
      TruePosC.threshold <- TruePosV.threshold <- TrueNegC.threshold <- TrueNegV.threshold <- FalsePosC.threshold <- FalsePosV.threshold <- FalseNegC.threshold <- FalseNegV.threshold <- list()
    }
    if("ER" %in% outputs){
      ErrorRateC.threshold <- ErrorRateV.threshold <- list() 
      ErrorRateCglobal.threshold <- ErrorRateVglobal.threshold <- list() 
    }
  }
  
  if((nNoBin==0) & ("AUC" %in% outputs)){
    AUCc <- list() 
    AUCv <- list() 
    
    AUCc.global <- list()
    AUCv.global <- list()
  }

      
  
  
  # step 4. iterations on permutated Y
  for(j in 1:(npermut+1)){
   
    set.seed(seed=j)
    
    # P.1. Permutations 
    nObsP <- NA
    if(j==1) {nObsP <- 0}
    if(j>1) {
      if(is.null(nbObsPermut)==TRUE){
        nObsP <- sample(x = nr, size = 1, replace=FALSE)
      }else{
        nObsP <- nbObsPermut
      }
    }
    
    Ypermut <- Y
    
    if((nObsP>0) == TRUE){
      for(o in 1:nObsP){
        indObsPermut                  <- sample(x = nr, size = 2, replace=FALSE)
        Y1                            <- Ypermut$tab[indObsPermut[1],]
        Y2                            <- Ypermut$tab[indObsPermut[2],]
        Ypermut$tab[indObsPermut[1],] <- Y2
        Ypermut$tab[indObsPermut[2],] <- Y1
      }
    }

    rownames(Ypermut$tab)        <- rownames(Y$tab)
    
    # P.2. Correlations/similarities initial and permutated Y block
    res$cor.YYpermut.values[j,2:(q+1)]    <- sapply(1:q, function(Q) cor(Y$tab[,Q], Ypermut$tab[,Q]))
    res$prct.Ychange.values[j,2:(q+1)]    <- sapply(1:q, function(Q) (sum(Y$tab[,Q] != Ypermut$tab[,Q])/nr)) 
    res$RV.YYpermut.values[j,2]           <- coeffRV(Y$tab,Ypermut$tab)$rv
    res$prctGlob.Ychange.values[j,2]      <- sum((Y$tab - Ypermut$tab)^2)/(nr*q)
 
    # P.3. Preparation of the parallelized processing
    #nodes <- detectCores()
    cl    <- makeCluster(cpus, type="PSOCK") # initialisation 
    registerDoParallel(cl) #  cl scripts creation
    on.exit(stopCluster(cl))
  
    # P.4. start of the parallelized iterations (repetitions / cross validation)
    resForeach <- foreach(i = 1:nrepet, .export=c ("mbplsda", "inertie", "ginv"), .packages=c("ade4","pROC"), .errorhandling="remove") %dopar%{
 
      set.seed(seed=i)
      
      ## I.1. Dividing X and Y into calibration (Xc, Yc) and validation (Xv, Yv) datasets. 
      s  <- sample(x = nr, size = Nc)  
      Xc <- X[, s, ]
      Xv <- X[, -s, ]
      Yc <- Ypermut[s, ]
      Yv <- Ypermut[-s, ]
      
      rnamesXc <- row.names(Xc)
      rnamesXv <- row.names(Xv)
      rnamesYc <- row.names(Yc$tab)
      rnamesYv <- row.names(Yv$tab)
      
      ## number of "1" by Y categorie
      nbY1c <- sapply(1:q, function(g) sum(Yc$tab[, g] == 1))
      nbY1v <- sapply(1:q, function(g) sum(Yv$tab[, g] == 1))
      
      
      ## I.2. Application of mbplsda on calibration/validation datasets  
      rescal <- do.call(method, list(dudiY = Yc, ktabX = Xc, scale = scale, option = option, scannf = FALSE, nf = optdim))# nf=h remplace par nf=optdim
      resval <- do.call(method, list(dudiY = Yv, ktabX = Xv, scale = scale, option = option, scannf = FALSE, nf = optdim)) # nf=h remplace par nf=optdim                  
      
    
      ## I.3. variable matrix
      ## Xc brute 
      Xc.mat <- cbind.data.frame(unclass(Xc)[1:nblo], stringsAsFactors = TRUE)
      
      ## Means and biased sd of Xc variables (to use with mbpls in ade4)
      rescal$meanX        <- colMeans(Xc.mat) # equivalent to rescal$meanX
      rescal$sdX          <- apply(Xc.mat, 2, sd) * sqrt((Nc-1)/Nc)  # biaised sd, equivalent to rescal$sdX
      
      ## Xv: Xv brute, Xv centred considering Xc means, Xv centred reduced weighted considering means, sd and intertia of Xc blocks
      Xv.raw <- cbind.data.frame(unclass(Xv)[1:nblo], stringsAsFactors = TRUE)
      
      Xv.c   <- sweep(Xv.raw, 2, rescal$meanX, FUN="-") # equivalent to Xv.raw - rep(rescal$meanX, each=Nv)
      
      if(scale==TRUE){
        Xv.cr  <- sweep(Xv.c, 2, rescal$sdX, FUN="/") # equivalent to Xv.c / rep(rescal$sdX, each=Nv)
        if(option=="uniform"){
          Xv.crw <- Xv.cr * sqrt(matrix(rep(rescal$X.cw, each=Nv), nrow=Nv)) # rescal$X.cw = column weights
        }else{
          Xv.crw <- Xv.cr
        }
      }
      
      if(scale==FALSE){
        Xv.cr  <- Xv.c
        if(option=="uniform"){
          Xv.crw <- Xv.cr * sqrt(matrix(rep(rescal$X.cw, each=Nv), nrow=Nv)) # rescal$X.cw = column weights
        }else{
          Xv.crw <- Xv.cr
        }
      }

      # I.4. preparation of outputs obtained with optdim
      ## index of predicted categories
      classYc           <- matrix(NA,Nc,ncol=3*Ky,dimnames=list(rnamesYc, c(paste0("Ymax_Var",1:Ky),paste0("Ycentroid_Var",1:Ky), paste0("Ythreshold_Var",1:Ky))))
      classYv           <- matrix(NA,Nv,ncol=3*Ky,dimnames=list(rnamesYv, c(paste0("Ymax_Var",1:Ky),paste0("Ycentroid_Var",1:Ky), paste0("Ythreshold_Var",1:Ky))))

      if("max" %in% algo){
        ## disjonctive table of predicted modalities with max
        ClassPredYc.max <- matrix(0,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
        ClassPredYv.max <- matrix(0,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
      }
      
      if("gravity" %in% algo){
        ## disjonctive table of predicted modalities with gravity
        ClassPredYc.gravity <- matrix(0,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
        ClassPredYv.gravity <- matrix(0,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
      }
      
      if("threshold" %in% algo){
        ## disjonctive table of predicted modalities with threshold
        ClassPredYc.threshold <- matrix(NA,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
        ClassPredYv.threshold <- matrix(NA,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
      }
      
      
      ## error rates by variable and overall
  	  RMSE_ErrorRates <- list(NULL)
  	  
  	  if("ER" %in% outputs){
  	    ErrorRateGlobal <- matrix(NA,nrow=6,ncol=(Ky+1),
  	                              dimnames=list(c("ErrRateCmaxGlobal","ErrRateVmaxGlobal","ErrRateCgravityGlobal","ErrRateVgravityGlobal","ErrRateCthresholdGlobal","ErrRateVthresholdGlobal"), c(paste0("Var",1:Ky),"global")))
  	  }
  	  if((nNoBin==0) & ("AUC" %in% outputs)){
  	    ## AUC by variable and overall
        AUCglobal       <- matrix(NA,nrow=2,ncol=(Ky+1), dimnames=list(c("AUCcGlobal","AUCvGlobal"), c(paste0("Var",1:Ky),"Mean")))
  	  }
  	  
      
      # I.5. coefficients appliable on raw Xc 
      XYcoef.raw.cal    <- sapply(rescal$XYcoef.raw, function(x) x[, optdim])
      intercept.raw.cal <- sapply(rescal$intercept.raw, function(x) x[, optdim])
      

      # I.6. Yc and Yv predicted by optdim
      predYc        <- matrix(rep(intercept.raw.cal, each=Nc), ncol=q) + as.matrix(Xc.mat) %*% XYcoef.raw.cal
      predYv        <- matrix(rep(intercept.raw.cal, each=Nv), ncol=q) + as.matrix(Xv.raw) %*% XYcoef.raw.cal
      colnames(predYc) <- colnames(predYv) <- cnames

      
      # I.7. RMSE
  	  if("RMSE" %in% outputs){
  	    residYc                       <- as.matrix(Yc$tab) - predYc
  	    RMSE_ErrorRates$RMSECglobal   <- sqrt(sum(residYc^2) / (Nc * q))
  	    RMSE_ErrorRates$RMSEC         <- sqrt(apply((residYc^2), mean, MARGIN = 2, na.rm=TRUE))
  	    residYv                       <- as.matrix(Yv$tab) - predYv
  	    RMSE_ErrorRates$RMSEVglobal   <- sqrt(sum(residYv^2) / (Nv * q))
  	    RMSE_ErrorRates$RMSEV         <- sqrt(apply((residYv^2), mean, MARGIN = 2, na.rm=TRUE))
  	  }
      
      
      # I.8. predictions and errors with max
      if("max" %in% algo){
        for (k in 1:Ky){
          
          ## index of predicted modalities with max
          classYc[,k] <- sapply(1:Nc, function(n) which.max(predYc[n,Var == k]))
          classYv[,k] <- sapply(1:Nv, function(n) which.max(predYv[n,Var == k]))
          
          ## indicators of predicted modalities with max in a disjunctive table
          for(n in 1:Nc){
            if(k==1) (ClassPredYc.max[n,classYc[n,1]] <- 1)
            if(k>1)  (ClassPredYc.max[n,(sum(bloY[1:(k-1)])+classYc[n,k])] <- 1)
          }
          for(n in 1:Nv){
            if(k==1) (ClassPredYv.max[n,classYv[n,1]] <- 1)
            if(k>1)  (ClassPredYv.max[n,(sum(bloY[1:(k-1)])+classYv[n,k])] <- 1)
          }
        }
        
        ## confusion matrix and error with max
        if(("ConfMat" %in% outputs) | ("ER" %in% outputs)){
          RMSE_ErrorRates$TPcM  <- sapply(1:q, function(l)length(which(ClassPredYc.max[,l]==1 & Yc$tab[,l]==1))/Nc)
          RMSE_ErrorRates$TPvM  <- sapply(1:q, function(l)length(which(ClassPredYv.max[,l]==1 & Yv$tab[,l]==1))/Nv)
        }
        if("ConfMat" %in% outputs){
          RMSE_ErrorRates$TNcM  <- sapply(1:q, function(l)length(which(ClassPredYc.max[,l]==0 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FPcM  <- sapply(1:q, function(l)length(which(ClassPredYc.max[,l]==1 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FNcM  <- sapply(1:q, function(l)length(which(ClassPredYc.max[,l]==0 & Yc$tab[,l]==1))/Nc)
          
          RMSE_ErrorRates$TNvM  <- sapply(1:q, function(l)length(which(ClassPredYv.max[,l]==0 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FPvM  <- sapply(1:q, function(l)length(which(ClassPredYv.max[,l]==1 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FNvM  <- sapply(1:q, function(l)length(which(ClassPredYv.max[,l]==0 & Yv$tab[,l]==1))/Nv)
        }
        if("ER" %in% outputs){
          RMSE_ErrorRates$ERcM  <- sapply(1:q, function(l)(length(which(ClassPredYc.max[,l]==1 & Yc$tab[,l]==0)) + length(which(ClassPredYc.max[,l]==0 & Yc$tab[,l]==1)))/Nc)
          RMSE_ErrorRates$ERvM  <- sapply(1:q, function(l)(length(which(ClassPredYv.max[,l]==1 & Yv$tab[,l]==0)) + length(which(ClassPredYv.max[,l]==0 & Yv$tab[,l]==1)))/Nv)
        }
      }
      
      
      # I.9. predictions and errors with gravity
      if("gravity" %in% algo){
        ## barycenters of groups on global components
        
        Gravity <- matrix(0, nrow=q, ncol=optdim, dimnames=list(cnames, 1:optdim))
        
        if (optdim == 1){
          Gravity[,1] <- (sapply(1:q, function(g){
            if(nbY1c[g]>1) (mean(rescal$lX[Yc$tab[, g] == 1, 1]))
            else if(nbY1c[g]==1) (rescal$lX[Yc$tab[, g] == 1, 1])
            else if(nbY1c[g]==0) (NA)
          }))
        }else{
          Gravity[,1:optdim] <- (t(sapply(1:q, function(g) {
            if(nbY1c[g]>1) (apply(rescal$lX[Yc$tab[, g] == 1,1:optdim], 2, mean))
            else if(nbY1c[g]==1) (rescal$lX[Yc$tab[, g] == 1,1:optdim])
            else if(nbY1c[g]==0) (rep(NA,times=optdim))
          })))
        }
        
        ## distances to barycenters
        dist.eucl.gravity.Yc    <- sapply(1:q, function(g) apply((rescal$lX[,1:optdim] -  matrix(rep(Gravity[g, ], each = Nc), nrow = Nc))^2, 1, sum))
        dist.eucl.gravity.Yv    <- sapply(1:q, function(g) apply(((as.matrix(Xv.crw)%*%rescal$faX[,1:optdim]) -  matrix(rep(Gravity[g, ], each = Nv), nrow = Nv))^2, 1, sum))
        
        for (k in 1:Ky){
          
          ## index of predicted modalities with gravity
          classYc[,(k+Ky)] <- sapply(1:Nc, function(n) which.min(dist.eucl.gravity.Yc[n,Var == k]))
          classYv[,(k+Ky)] <- sapply(1:Nv, function(n) which.min(dist.eucl.gravity.Yv[n,Var == k]))
          
          ## indicators of predicted modalities with gravity in a disjunctive table
          for(n in 1:Nc){
            if(k==1) (ClassPredYc.gravity[n,classYc[n,(1+Ky)]] <- 1)
            if(k>1)  (ClassPredYc.gravity[n,(sum(bloY[1:(k-1)])+classYc[n,(k+Ky)])] <- 1)
          }
          for(n in 1:Nv){
            if(k==1) (ClassPredYv.gravity[n,classYv[n,(1+Ky)]] <- 1)
            if(k>1)  (ClassPredYv.gravity[n,(sum(bloY[1:(k-1)])+classYv[n,(k+Ky)])] <- 1)
          }
        }
        
        ## confusion matrix and error with gravity
        if(("ConfMat" %in% outputs) | ("ER" %in% outputs)){
          RMSE_ErrorRates$TPcG <- sapply(1:q, function(l)length(which(ClassPredYc.gravity[,l]==1 & Yc$tab[,l]==1))/Nc)
          RMSE_ErrorRates$TPvG <- sapply(1:q, function(l)length(which(ClassPredYv.gravity[,l]==1 & Yv$tab[,l]==1))/Nv)
        }
        if("ConfMat" %in% outputs){
          RMSE_ErrorRates$TNcG <- sapply(1:q, function(l)length(which(ClassPredYc.gravity[,l]==0 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FPcG <- sapply(1:q, function(l)length(which(ClassPredYc.gravity[,l]==1 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FNcG <- sapply(1:q, function(l)length(which(ClassPredYc.gravity[,l]==0 & Yc$tab[,l]==1))/Nc)
          
          RMSE_ErrorRates$TNvG <- sapply(1:q, function(l)length(which(ClassPredYv.gravity[,l]==0 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FPvG <- sapply(1:q, function(l)length(which(ClassPredYv.gravity[,l]==1 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FNvG <- sapply(1:q, function(l)length(which(ClassPredYv.gravity[,l]==0 & Yv$tab[,l]==1))/Nv)
        }
        if("ER" %in% outputs){
          RMSE_ErrorRates$ERcG <- sapply(1:q, function(l)(length(which(ClassPredYc.gravity[,l]==1 & Yc$tab[,l]==0)) + length(which(ClassPredYc.gravity[,l]==0 & Yc$tab[,l]==1)))/Nc)
          RMSE_ErrorRates$ERvG <- sapply(1:q, function(l)(length(which(ClassPredYv.gravity[,l]==1 & Yv$tab[,l]==0)) + length(which(ClassPredYv.gravity[,l]==0 & Yv$tab[,l]==1)))/Nv)
        }
      }
      
      
      # I.10. predictions and errors with threshold
      if("threshold" %in% algo){
        ## indicators of predicted modalities with threshold in a disjunctive table
        ClassPredYc.threshold[predYc>=threshold] <- 1
        ClassPredYv.threshold[predYv>=threshold] <- 1
        ClassPredYc.threshold[predYc<threshold]  <- 0
        ClassPredYv.threshold[predYv<threshold]  <- 0
        
        ## confusion matrix and error with threshold
        if(("ConfMat" %in% outputs) | ("ER" %in% outputs)){
          RMSE_ErrorRates$TPcT <- sapply(1:q, function(l)length(which(ClassPredYc.threshold[,l]==1 & Yc$tab[,l]==1))/Nc)
          RMSE_ErrorRates$TPvT <- sapply(1:q, function(l)length(which(ClassPredYv.threshold[,l]==1 & Yv$tab[,l]==1))/Nv)
        }
        if("ConfMat" %in% outputs){
          RMSE_ErrorRates$TNcT <- sapply(1:q, function(l)length(which(ClassPredYc.threshold[,l]==0 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FPcT <- sapply(1:q, function(l)length(which(ClassPredYc.threshold[,l]==1 & Yc$tab[,l]==0))/Nc)
          RMSE_ErrorRates$FNcT <- sapply(1:q, function(l)length(which(ClassPredYc.threshold[,l]==0 & Yc$tab[,l]==1))/Nc)
          
          RMSE_ErrorRates$TNvT <- sapply(1:q, function(l)length(which(ClassPredYv.threshold[,l]==0 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FPvT <- sapply(1:q, function(l)length(which(ClassPredYv.threshold[,l]==1 & Yv$tab[,l]==0))/Nv)
          RMSE_ErrorRates$FNvT <- sapply(1:q, function(l)length(which(ClassPredYv.threshold[,l]==0 & Yv$tab[,l]==1))/Nv)
        }
        if("ER" %in% outputs){
          RMSE_ErrorRates$ERcT <- sapply(1:q, function(l)(length(which(ClassPredYc.threshold[,l]==1 & Yc$tab[,l]==0)) + length(which(ClassPredYc.threshold[,l]==0 & Yc$tab[,l]==1)))/Nc)
          RMSE_ErrorRates$ERvT <- sapply(1:q, function(l)(length(which(ClassPredYv.threshold[,l]==1 & Yv$tab[,l]==0)) + length(which(ClassPredYv.threshold[,l]==0 & Yv$tab[,l]==1)))/Nv)
        }
        
        for (k in 1:Ky){
          ## index of predicted modalities with threshold
          classYc[,(2*Ky+k)] <- sapply(1:Nc, function(n) if(sum(ClassPredYc.threshold[n,Var == k], na.rm=T)==1) (which.max(ClassPredYc.threshold[n,Var == k])))
          classYv[,(2*Ky+k)] <- sapply(1:Nv, function(n) if(sum(ClassPredYc.threshold[n,Var == k], na.rm=T)==1) (which.max(ClassPredYv.threshold[n,Var == k])))      
        }
      }
      
      if((nNoBin==0) & ("AUC" %in% outputs)){
        ## AUC if the 2 categories of each Y variable are always in Yc and Yv
        if(sum(sum(nbY1v==0), sum(nbY1c==0))==0){
          RMSE_ErrorRates$AUCc <- sapply(1:q, function(l) auc(Yc$tab[,l],predYc[,l]))
          RMSE_ErrorRates$AUCv <- sapply(1:q, function(l) auc(Yv$tab[,l],predYv[,l]))
        }
      }
      
      ## I.11. global error rates
      
      if((nNoBin==0) & ("AUC" %in% outputs)){
        AUCglobal["AUCcGlobal",1:Ky]       <- RMSE_ErrorRates$AUCc[seq(from=1, to=q, by=2)]
        AUCglobal["AUCvGlobal",1:Ky]       <- RMSE_ErrorRates$AUCv[seq(from=1, to=q, by=2)]
        AUCglobal["AUCcGlobal","Mean"]     <- mean(RMSE_ErrorRates$AUCc[seq(from=1, to=q, by=2)],na.rm=TRUE)
        AUCglobal["AUCvGlobal","Mean"]     <- mean(RMSE_ErrorRates$AUCv[seq(from=1, to=q, by=2)],na.rm=TRUE)
        RMSE_ErrorRates$AUCcGlobal <- AUCglobal["AUCcGlobal",]
        RMSE_ErrorRates$AUCvGlobal <- AUCglobal["AUCvGlobal",]
      }
      
      if("ER" %in% outputs){ 
        if("max" %in% algo){
          ErrorRateGlobal["ErrRateCmaxGlobal",1:Ky]       <- sapply(1:Ky, function(k)(1-sum(RMSE_ErrorRates$TPcM[Var == k])))
          ErrorRateGlobal["ErrRateVmaxGlobal",1:Ky]       <- sapply(1:Ky, function(k)(1-sum(RMSE_ErrorRates$TPvM[Var == k])))
          ErrorRateGlobal["ErrRateCmaxGlobal","global"]   <- 1-sum((apply(((ClassPredYc.max - Yc$tab)^2), sum, MARGIN = 1))==0)/Nc
          ErrorRateGlobal["ErrRateVmaxGlobal","global"]   <- 1-sum((apply(((ClassPredYv.max - Yv$tab)^2), sum, MARGIN = 1))==0)/Nv
          RMSE_ErrorRates$ERcMglobal <- ErrorRateGlobal["ErrRateCmaxGlobal",]
          RMSE_ErrorRates$ERvMglobal <- ErrorRateGlobal["ErrRateVmaxGlobal",]
        }
        
        if("gravity" %in% algo){
          ErrorRateGlobal["ErrRateCgravityGlobal",1:Ky]     <- sapply(1:Ky, function(k)(1-sum(RMSE_ErrorRates$TPcG[Var == k])))
          ErrorRateGlobal["ErrRateVgravityGlobal",1:Ky]     <- sapply(1:Ky, function(k)(1-sum(RMSE_ErrorRates$TPvG[Var == k])))
          ErrorRateGlobal["ErrRateCgravityGlobal","global"] <- 1-sum((apply(((ClassPredYc.gravity - Yc$tab)^2), sum, MARGIN = 1))==0)/Nc
          ErrorRateGlobal["ErrRateVgravityGlobal","global"] <- 1-sum((apply(((ClassPredYv.gravity - Yv$tab)^2), sum, MARGIN = 1))==0)/Nv
          RMSE_ErrorRates$ERcGglobal <- ErrorRateGlobal["ErrRateCgravityGlobal",]
          RMSE_ErrorRates$ERvGglobal <- ErrorRateGlobal["ErrRateVgravityGlobal",]
        }
        
        if("threshold" %in% algo){
          ErrorRateGlobal["ErrRateCthresholdGlobal",1:Ky]     <- sapply(1:Ky, function(k)(1-sum((apply(((ClassPredYc.threshold[,Var == k] - Yc$tab[,Var == k])^2), sum, MARGIN = 1))==0)/Nc))
          ErrorRateGlobal["ErrRateVthresholdGlobal",1:Ky]     <- sapply(1:Ky, function(k)(1-sum((apply(((ClassPredYv.threshold[,Var == k] - Yv$tab[,Var == k])^2), sum, MARGIN = 1))==0)/Nv))
          ErrorRateGlobal["ErrRateCthresholdGlobal","global"] <- 1-sum((apply(((ClassPredYc.threshold - Yc$tab)^2), sum, MARGIN = 1))==0)/Nc
          ErrorRateGlobal["ErrRateVthresholdGlobal","global"] <- 1-sum((apply(((ClassPredYv.threshold - Yv$tab)^2), sum, MARGIN = 1))==0)/Nv
          RMSE_ErrorRates$ERcTglobal <- ErrorRateGlobal["ErrRateCthresholdGlobal",]
          RMSE_ErrorRates$ERvTglobal <- ErrorRateGlobal["ErrRateVthresholdGlobal",]
        }
      }
      
 
      ## I.12. list of foreach results

      RMSE_ErrorRates
    }  
    
    # step 4. End of the parallelized loops
  
    stopCluster(cl)
    on.exit(stopCluster)
#    resForeach
    
    nrepetFE       <- length(resForeach)
    if((nrepetFE<1.5)|(is.null(nrepetFE)==TRUE)){
      stop("No adjustement of models")
    }
    #res$TRUEnrepet <- nrepetFE
    
    # step 5. prepare outputs 

    dimlabR <- paste("repet", (1 : nrepetFE), sep = "")
    
    ## concatenated results of cv repetitions on the j th permutated Y block
    if("max" %in% algo){
      if("ConfMat" %in% outputs){
        TPcMall <- TPvMall <- TNcMall <- TNvMall <- FPcMall <- FPvMall <- FNcMall <- FNvMall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
      }
      if("ER" %in% outputs){
        ERcMall <- ERvMall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
        ERcMglobalall <- ERvMglobalall <- matrix(NA, nrow = nrepetFE, ncol = (Ky+1), dimnames=list(dimlabR,c(paste0("Var",1:Ky),"global")))
      }
    }
    if("gravity" %in% algo){
      if("ConfMat" %in% outputs){
        TPcGall <- TPvGall <- TNcGall <- TNvGall <- FPcGall <- FPvGall <- FNcGall <- FNvGall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
      }
      if("ER" %in% outputs){
        ERcGall <- ERvGall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
        ERcGglobalall <- ERvGglobalall <- matrix(NA, nrow = nrepetFE, ncol = (Ky+1), dimnames=list(dimlabR,c(paste0("Var",1:Ky),"global")))
      }
    }
    
    if("threshold" %in% algo){
      if("ConfMat" %in% outputs){
        TPcTall <- TPvTall <- TNcTall <- TNvTall <- FPcTall <- FPvTall <- FNcTall <- FNvTall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
      }
      if("ER" %in% outputs){
        ERcTall <- ERvTall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
        ERcTglobalall <- ERvTglobalall <- matrix(NA, nrow = nrepetFE, ncol = (Ky+1), dimnames=list(dimlabR,c(paste0("Var",1:Ky),"global")))
      }
    }
    
    if((nNoBin == 0) & ("AUC" %in% outputs)){
      AUCcall <- AUCvall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
      AUCcglobalall <- AUCvglobalall <- matrix(NA, nrow = nrepetFE, ncol = (Ky+1), dimnames=list(dimlabR,c(paste0("Var",1:Ky),"Mean")))
    }
    if("RMSE" %in% outputs){
      RMSECglobalall <- RMSEVglobalall <- rep(NA, nrepetFE)
      RMSECall <- RMSEVall <- matrix(NA, nrow=nrepetFE, ncol=q, dimnames=list(dimlabR,cnames))
    }
    

    # step 6. means on repetitions
    ### concatenated results 

    for(i in 1:nrepetFE){
      if("RMSE" %in% outputs){
        RMSECglobalall[i] <- (resForeach[[i]][["RMSECglobal"]])
        RMSEVglobalall[i] <- (resForeach[[i]][["RMSEVglobal"]])
        RMSECall[i,] <- (resForeach[[i]][["RMSEC"]])
        RMSEVall[i,] <- (resForeach[[i]][["RMSEV"]])
      }
      
      if("max" %in% algo){
        if("ConfMat" %in% outputs){
          TPcMall[i,] <- (resForeach[[i]][["TPcM"]])
          TPvMall[i,] <- (resForeach[[i]][["TPvM"]])
          TNcMall[i,] <- (resForeach[[i]][["TNcM"]])
          TNvMall[i,] <- (resForeach[[i]][["TNvM"]])
          FPcMall[i,] <- (resForeach[[i]][["FPcM"]])
          FPvMall[i,] <- (resForeach[[i]][["FPvM"]])
          FNcMall[i,] <- (resForeach[[i]][["FNcM"]])
          FNvMall[i,] <- (resForeach[[i]][["FNvM"]])
        }
        if("ER" %in% outputs){
          ERcMall[i,] <- (resForeach[[i]][["ERcM"]])
          ERvMall[i,] <- (resForeach[[i]][["ERvM"]])
          ERcMglobalall[i,] <- resForeach[[i]][["ERcMglobal"]]
          ERvMglobalall[i,] <- resForeach[[i]][["ERvMglobal"]]
        }
      }
      
      if("gravity" %in% algo){
        if("ConfMat" %in% outputs){
          TPcGall[i,] <- (resForeach[[i]][["TPcG"]])
          TPvGall[i,] <- (resForeach[[i]][["TPvG"]])
          TNcGall[i,] <- (resForeach[[i]][["TNcG"]])
          TNvGall[i,] <- (resForeach[[i]][["TNvG"]])
          FPcGall[i,] <- (resForeach[[i]][["FPcG"]])
          FPvGall[i,] <- (resForeach[[i]][["FPvG"]])
          FNcGall[i,] <- (resForeach[[i]][["FNcG"]])
          FNvGall[i,] <- (resForeach[[i]][["FNvG"]])
        }
        if("ER" %in% outputs){
          ERcGall[i,] <- (resForeach[[i]][["ERcG"]])
          ERvGall[i,] <- (resForeach[[i]][["ERvG"]])
          ERcGglobalall[i,] <- resForeach[[i]][["ERcGglobal"]]
          ERvGglobalall[i,] <- resForeach[[i]][["ERvGglobal"]]
        }
      }
      
      if("threshold" %in% algo){
        if("ConfMat" %in% outputs){
          TPcTall[i,] <- (resForeach[[i]][["TPcT"]])
          TPvTall[i,] <- (resForeach[[i]][["TPvT"]])
          TNcTall[i,] <- (resForeach[[i]][["TNcT"]])
          TNvTall[i,] <- (resForeach[[i]][["TNvT"]])
          FPcTall[i,] <- (resForeach[[i]][["FPcT"]])
          FPvTall[i,] <- (resForeach[[i]][["FPvT"]])
          FNcTall[i,] <- (resForeach[[i]][["FNcT"]])
          FNvTall[i,] <- (resForeach[[i]][["FNvT"]])
        }
        if("ER" %in% outputs){
          ERcTall[i,] <- (resForeach[[i]][["ERcT"]])
          ERvTall[i,] <- (resForeach[[i]][["ERvT"]])
          ERcTglobalall[i,] <- resForeach[[i]][["ERcTglobal"]]
          ERvTglobalall[i,] <- resForeach[[i]][["ERvTglobal"]]
        }
      }
      
      if((nNoBin==0) & ("AUC" %in% outputs)){  
        AUCcall[i,] <- (resForeach[[i]][["AUCc"]])
        AUCvall[i,] <- (resForeach[[i]][["AUCv"]])
        AUCcglobalall[i,] <- resForeach[[i]][["AUCcGlobal"]]
        AUCvglobalall[i,] <- resForeach[[i]][["AUCvGlobal"]]
      }
    }
    
    

    ## functions for stat
    IC95 <- function(m){ # m is a vector
      rt <- c(rep(NA,2))
      if((sd(m, na.rm=TRUE)!=0) & (is.na(sd(m, na.rm=TRUE))==FALSE)){
        testt <- t.test(m, conf.level = 0.95)
        rt <- c(round(testt$conf.int[1],5),round(testt$conf.int[2],5))
      }
      return(rt)
    }
    
    
    stat.desc   <- function (x){ # x is a matrix, columns are variables, raws are repetitions
      nombre    <- colSums(!is.na(x))
      moy       <- round(colMeans(x,na.rm = TRUE),5)
      etype     <- round(apply(x,2,sd, na.rm=TRUE),5)
      quartiles <- round(t(apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),5)
      IC        <- t(apply(x, 2, IC95))
      result    <- cbind.data.frame(nombre, moy, etype, IC, quartiles, stringsAsFactors = TRUE)
      colnames(result) <- c("nb", "mean", "sd", "95CIinf", "95CIsup","Q2.5", "median", "Q97.5")
      rownames(result) <- colnames(x)
      return(result)
    }
    
    
    ### means on the nrepetFE repetitions  
    dataP <- rep(dimlabP[j],q)
    
    if("RMSE" %in% outputs){
      RMSEC.P[[j]]       <- cbind(cname=cnames,dataP,stat.desc(RMSECall))
      RMSECglobal.P[[j]] <- cbind(dimlabP[j],stat.desc(as.matrix(RMSECglobalall)))
      RMSEV.P[[j]]       <- cbind(cname=cnames,dataP,stat.desc(RMSEVall))
      RMSEVglobal.P[[j]] <- cbind(dimlabP[j],stat.desc(as.matrix(RMSEVglobalall)))
      colnames(RMSECglobal.P[[j]])[1] <- colnames(RMSEVglobal.P[[j]])[1] <- "dataP"
    }
    
    if("max" %in% algo){
      if("ConfMat" %in% outputs){
        TruePosC.max[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPcMall))
        TruePosV.max[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPvMall))
        TrueNegC.max[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNcMall))
        TrueNegV.max[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNvMall))
        FalsePosC.max[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPcMall))
        FalsePosV.max[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPvMall))
        FalseNegC.max[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNcMall))
        FalseNegV.max[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNvMall))
      }
      if("ER" %in% outputs){
        ErrorRateC.max[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERcMall))
        ErrorRateV.max[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERvMall))
        ErrorRateCglobal.max[[j]] <- cbind(colnames(ERcMglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERcMglobalall))
        ErrorRateVglobal.max[[j]] <- cbind(colnames(ERvMglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERvMglobalall))
        colnames(ErrorRateCglobal.max[[j]])[1:2] <- colnames(ErrorRateVglobal.max[[j]])[1:2] <- c("variable","dataP")
      }
    }
    
    if("gravity" %in% algo){
      if("ConfMat" %in% outputs){
        TruePosC.gravity[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPcGall))
        TruePosV.gravity[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPvGall))
        TrueNegC.gravity[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNcGall))
        TrueNegV.gravity[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNvGall))
        FalsePosC.gravity[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPcGall))
        FalsePosV.gravity[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPvGall))
        FalseNegC.gravity[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNcGall))
        FalseNegV.gravity[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNvGall))
      }
      if("ER" %in% outputs){
        ErrorRateC.gravity[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERcGall))
        ErrorRateV.gravity[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERvGall))
        ErrorRateCglobal.gravity[[j]] <- cbind(colnames(ERcGglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERcGglobalall))
        ErrorRateVglobal.gravity[[j]] <- cbind(colnames(ERvGglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERvGglobalall))
        colnames(ErrorRateCglobal.gravity[[j]])[1:2] <- colnames(ErrorRateVglobal.gravity[[j]])[1:2] <- c("variable","dataP")
      }
    }
    
    if("threshold" %in% algo){
      if("ConfMat" %in% outputs){
        TruePosC.threshold[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPcTall))
        TruePosV.threshold[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TPvTall))
        TrueNegC.threshold[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNcTall))
        TrueNegV.threshold[[j]]   <- cbind(cname=cnames,dataP,stat.desc(TNvTall))
        FalsePosC.threshold[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPcTall))
        FalsePosV.threshold[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FPvTall))
        FalseNegC.threshold[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNcTall))
        FalseNegV.threshold[[j]]  <- cbind(cname=cnames,dataP,stat.desc(FNvTall))
      }
      if("ER" %in% outputs){
        ErrorRateC.threshold[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERcTall))
        ErrorRateV.threshold[[j]] <- cbind(cname=cnames,dataP,stat.desc(ERvTall))
        ErrorRateCglobal.threshold[[j]] <- cbind(colnames(ERcTglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERcTglobalall))
        ErrorRateVglobal.threshold[[j]] <- cbind(colnames(ERvTglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(ERvTglobalall))
        colnames(ErrorRateCglobal.threshold[[j]])[1:2] <- colnames(ErrorRateVglobal.threshold[[j]])[1:2] <- c("variable","dataP")
      }
    }
    
    if((nNoBin==0) & ("AUC" %in% outputs)){
      AUCc[[j]] <- cbind(cname=cnames,dataP,stat.desc(AUCcall))
      AUCv[[j]] <- cbind(cname=cnames,dataP,stat.desc(AUCvall))
      AUCc.global[[j]] <- cbind(colnames(AUCcglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(AUCcglobalall))
      AUCv.global[[j]] <- cbind(colnames(AUCvglobalall), rep(dimlabP[j],(Ky+1)),stat.desc(AUCvglobalall))
      colnames(AUCc.global[[j]])[1:2] <- colnames(AUCv.global[[j]])[1:2] <- c("variable","dataP")
    }

    
  } # end of the loop on permutated Y blocks
 
  # means on permutated Y blocks
  res$descrYperm <- data.frame(descrYperm=c("RV.YYpermut",paste0("cor.YYpermut_",cnames), "prctGlob.YYpermut", paste0("prct.YYpermut_",cnames)), 
                               rbind(stat.desc(as.matrix(res$RV.YYpermut.values[2:(npermut+1),2])), stat.desc(as.matrix(res$cor.YYpermut.values[2:(npermut+1),2:(q+1)])),
                                     stat.desc(as.matrix(res$prctGlob.Ychange.values[2:(npermut+1),2])), stat.desc(as.matrix(res$prct.Ychange.values[2:(npermut+1),2:(q+1)]))))
  rownames(res$descrYperm) <- NULL

  # concatenation of results on permutated Y blocks
  
  ordreDF <- function(DFrame, colonne){
    DFrame <- DFrame[order(DFrame[,colonne]),]
    DFrame
  }
  if("RMSE" %in% outputs){
    res$RMSEc         <- do.call("rbind", RMSEC.P)
    res$RMSEc         <- ordreDF(res$RMSEc,"cname")
    res$RMSEc.global  <- do.call("rbind", RMSECglobal.P)
    res$RMSEv         <- do.call("rbind", RMSEV.P )
    res$RMSEv         <- ordreDF(res$RMSEv,"cname")
    res$RMSEVc.global <- do.call("rbind", RMSEVglobal.P)
    rownames(res$RMSEc) <- rownames(res$RMSEc.global) <- rownames(res$RMSEv) <- rownames(res$RMSEv.global) <- NULL
  }
  
  if("max" %in% algo){
    if("ConfMat" %in% outputs){
      res$TruePosC.max   <- do.call("rbind", TruePosC.max)
      res$TruePosC.max   <- ordreDF(res$TruePosC.max,"cname")
      res$TruePosV.max   <- do.call("rbind", TruePosV.max)
      res$TruePosV.max   <- ordreDF(res$TruePosV.max,"cname")
      res$TrueNegC.max   <- do.call("rbind", TrueNegC.max)
      res$TrueNegC.max   <- ordreDF(res$TrueNegC.max,"cname")
      res$TrueNegV.max   <- do.call("rbind", TrueNegV.max)
      res$TrueNegV.max   <- ordreDF(res$TrueNegV.max,"cname")
      res$FalsePosC.max  <- do.call("rbind", FalsePosC.max)
      res$FalsePosC.max  <- ordreDF(res$FalsePosC.max,"cname")
      res$FalsePosV.max  <- do.call("rbind", FalsePosV.max)
      res$FalsePosV.max  <- ordreDF(res$FalsePosV.max,"cname")
      res$FalseNegC.max  <- do.call("rbind", FalseNegC.max)
      res$FalseNegC.max  <- ordreDF(res$FalseNegC.max,"cname")
      res$FalseNegV.max  <- do.call("rbind", FalseNegV.max)
      res$FalseNegV.max  <- ordreDF(res$FalseNegV.max,"cname")
      
      rownames(res$TruePosC.max) <- rownames(res$TrueNegC.max) <- rownames(res$FalsePosC.max) <- rownames(res$FalseNegC.max) <- NULL
      rownames(res$TruePosV.max) <- rownames(res$TrueNegV.max) <- rownames(res$FalsePosV.max) <- rownames(res$FalseNegV.max) <- NULL
    }
    if("ER" %in% outputs){
      res$ErrorRateC.max <- do.call("rbind", ErrorRateC.max)
      res$ErrorRateC.max <- ordreDF(res$ErrorRateC.max,"cname")
      res$ErrorRateV.max <- do.call("rbind", ErrorRateV.max)
      res$ErrorRateV.max <- ordreDF(res$ErrorRateV.max,"cname")
      res$ErrorRateCglobal.max <- do.call("rbind", ErrorRateCglobal.max)
      res$ErrorRateCglobal.max <- ordreDF(res$ErrorRateCglobal.max,"variable")
      res$ErrorRateVglobal.max <- do.call("rbind", ErrorRateVglobal.max)
      res$ErrorRateVglobal.max <- ordreDF(res$ErrorRateVglobal.max,"variable")
      
      rownames(res$ErrorRateC.max) <- rownames(res$ErrorRateCglobal.max) <- NULL
      rownames(res$ErrorRateV.max) <- rownames(res$ErrorRateVglobal.max) <- NULL
    }
  }
  
  if("gravity" %in% algo){
    if("ConfMat" %in% outputs){
      res$TruePosC.gravity   <- do.call("rbind", TruePosC.gravity)
      res$TruePosC.gravity   <- ordreDF(res$TruePosC.gravity,"cname")
      res$TruePosV.gravity   <- do.call("rbind", TruePosV.gravity)
      res$TruePosV.gravity   <- ordreDF(res$TruePosV.gravity,"cname")
      res$TrueNegC.gravity   <- do.call("rbind", TrueNegC.gravity)
      res$TrueNegC.gravity   <- ordreDF(res$TrueNegC.gravity,"cname")
      res$TrueNegV.gravity   <- do.call("rbind", TrueNegV.gravity)
      res$TrueNegV.gravity   <- ordreDF(res$TrueNegV.gravity,"cname")
      res$FalsePosC.gravity  <- do.call("rbind", FalsePosC.gravity)
      res$FalsePosC.gravity  <- ordreDF(res$FalsePosC.gravity,"cname")
      res$FalsePosV.gravity  <- do.call("rbind", FalsePosV.gravity)
      res$FalsePosV.gravity  <- ordreDF(res$FalsePosV.gravity,"cname")
      res$FalseNegC.gravity  <- do.call("rbind", FalseNegC.gravity)
      res$FalseNegC.gravity  <- ordreDF(res$FalseNegC.gravity,"cname")
      res$FalseNegV.gravity  <- do.call("rbind", FalseNegV.gravity)
      res$FalseNegV.gravity  <- ordreDF(res$FalseNegV.gravity,"cname")
      
      rownames(res$TruePosC.gravity) <- rownames(res$TrueNegC.gravity) <- rownames(res$FalsePosC.gravity) <- rownames(res$FalseNegC.gravity) <- NULL
      rownames(res$TruePosV.gravity) <- rownames(res$TrueNegV.gravity) <- rownames(res$FalsePosV.gravity) <- rownames(res$FalseNegV.gravity) <- NULL
    }
    if("ER" %in% outputs){
      res$ErrorRateC.gravity <- do.call("rbind", ErrorRateC.gravity)
      res$ErrorRateC.gravity <- ordreDF(res$ErrorRateC.gravity,"cname")
      res$ErrorRateV.gravity <- do.call("rbind", ErrorRateV.gravity)
      res$ErrorRateV.gravity <- ordreDF(res$ErrorRateV.gravity,"cname")
      res$ErrorRateCglobal.gravity <- do.call("rbind", ErrorRateCglobal.gravity)
      res$ErrorRateCglobal.gravity <- ordreDF(res$ErrorRateCglobal.gravity,"variable")
      res$ErrorRateVglobal.gravity <- do.call("rbind", ErrorRateVglobal.gravity)
      res$ErrorRateVglobal.gravity <- ordreDF(res$ErrorRateVglobal.gravity,"variable")
      
      rownames(res$ErrorRateC.gravity) <- rownames(res$ErrorRateCglobal.gravity) <- NULL
      rownames(res$ErrorRateV.gravity) <- rownames(res$ErrorRateVglobal.gravity) <- NULL
    }
  }
  
  if("threshold" %in% algo){
    if("ConfMat" %in% outputs){
      res$TruePosC.threshold   <- do.call("rbind", TruePosC.threshold)
      res$TruePosC.threshold   <- ordreDF(res$TruePosC.threshold,"cname")
      res$TruePosV.threshold   <- do.call("rbind", TruePosV.threshold)
      res$TruePosV.threshold   <- ordreDF(res$TruePosV.threshold,"cname")
      res$TrueNegC.threshold   <- do.call("rbind", TrueNegC.threshold)
      res$TrueNegC.threshold   <- ordreDF(res$TrueNegC.threshold,"cname")
      res$TrueNegV.threshold   <- do.call("rbind", TrueNegV.threshold)
      res$TrueNegV.threshold   <- ordreDF(res$TrueNegV.threshold,"cname")
      res$FalsePosC.threshold  <- do.call("rbind", FalsePosC.threshold)
      res$FalsePosC.threshold  <- ordreDF(res$FalsePosC.threshold,"cname")
      res$FalsePosV.threshold  <- do.call("rbind", FalsePosV.threshold)
      res$FalsePosV.threshold  <- ordreDF(res$FalsePosV.threshold,"cname")
      res$FalseNegC.threshold  <- do.call("rbind", FalseNegC.threshold)
      res$FalseNegC.threshold  <- ordreDF(res$FalseNegC.threshold,"cname")
      res$FalseNegV.threshold  <- do.call("rbind", FalseNegV.threshold)
      res$FalseNegV.threshold  <- ordreDF(res$FalseNegV.threshold,"cname")
      res$ErrorRateC.threshold <- do.call("rbind", ErrorRateC.threshold)
      res$ErrorRateC.threshold <- ordreDF(res$ErrorRateC.threshold,"cname")
      res$ErrorRateV.threshold <- do.call("rbind", ErrorRateV.threshold)
      res$ErrorRateV.threshold <- ordreDF(res$ErrorRateV.threshold,"cname")
      
      rownames(res$TruePosC.threshold) <- rownames(res$TrueNegC.threshold) <- rownames(res$FalsePosC.threshold) <- rownames(res$FalseNegC.threshold) <- NULL
      rownames(res$TruePosV.threshold) <- rownames(res$TrueNegV.threshold) <- rownames(res$FalsePosV.threshold) <- rownames(res$FalseNegV.threshold) <- NULL
    }
    if("ER" %in% outputs){
      res$ErrorRateCglobal.threshold <- do.call("rbind", ErrorRateCglobal.threshold)
      res$ErrorRateCglobal.threshold <- ordreDF(res$ErrorRateCglobal.threshold,"variable")
      res$ErrorRateVglobal.threshold <- do.call("rbind", ErrorRateVglobal.threshold)
      res$ErrorRateVglobal.threshold <- ordreDF(res$ErrorRateVglobal.threshold,"variable")
      
      rownames(res$ErrorRateC.threshold) <- rownames(res$ErrorRateCglobal.threshold) <- NULL
      rownames(res$ErrorRateV.threshold) <- rownames(res$ErrorRateVglobal.threshold) <- NULL
    }
  }
  
  if((nNoBin==0) & ("AUC" %in% outputs)){
    res$AUCc <- do.call("rbind", AUCc)
    res$AUCc <- ordreDF(res$AUCc,"cname")
    res$AUCv <- do.call("rbind", AUCv)
    res$AUCv <- ordreDF(res$AUCv,"cname")
    res$AUCc.global <- do.call("rbind", AUCc.global)
    res$AUCc.global <- ordreDF(res$AUCc.global,"variable")
    res$AUCv.global <- do.call("rbind", AUCv.global)
    res$AUCv.global <- ordreDF(res$AUCv.global,"variable")
    
    rownames(res$AUCc) <- rownames(res$AUCv) <- rownames(res$AUCc.global) <- rownames(res$AUCv.global) <- NULL
  }
  
  
  # regressions 
  if(("ER" %in% outputs) | ("AUC" %in% outputs)){
    res$reg.GlobalRes_prctYchange <- matrix(NA, nrow=4, ncol=2, dimnames = list(c("ERvGlobal.max","ERvGlobal.gravity","ERvGlobal.threshold","AUCvGlobal"),c("coeff","pvalue")))
  }
  if("ER" %in% outputs){
    if("max" %in% algo){
      regressionM <- lm(res$ErrorRateVglobal.max[(res$ErrorRateVglobal.max$variable=="global" & res$ErrorRateVglobal.max$dataP %in% (dimlabP[2:(npermut+1)])),"mean"] ~ res$prctGlob.Ychange.values[2:(npermut+1),2])
      res$reg.GlobalRes_prctYchange["ERvGlobal.max","coeff"]  <- regressionM$coefficients[2]
      res$reg.GlobalRes_prctYchange["ERvGlobal.max","pvalue"] <- anova(regressionM)[1,"Pr(>F)"]
    }
    
  	if("gravity" %in% algo){
  	  regressionG <- lm(res$ErrorRateVglobal.gravity[(res$ErrorRateVglobal.gravity$variable=="global" & res$ErrorRateVglobal.gravity$dataP %in% (dimlabP[2:(npermut+1)])),"mean"] ~ res$prctGlob.Ychange.values[2:(npermut+1),2])
  	  res$reg.GlobalRes_prctYchange["ERvGlobal.gravity","coeff"]  <- regressionG$coefficients[2]
  	  res$reg.GlobalRes_prctYchange["ERvGlobal.gravity","pvalue"] <- anova(regressionG)[1,"Pr(>F)"]
  	}
    
    if("threshold" %in% algo){
      regressionT <- lm(res$ErrorRateVglobal.threshold[(res$ErrorRateVglobal.threshold$variable=="global" & res$ErrorRateVglobal.threshold$dataP %in% (dimlabP[2:(npermut+1)])),"mean"]	~ res$prctGlob.Ychange.values[2:(npermut+1),2])
      res$reg.GlobalRes_prctYchange["ERvGlobal.threshold","coeff"]  <- regressionT$coefficients[2]
      res$reg.GlobalRes_prctYchange["ERvGlobal.threshold","pvalue"] <- anova(regressionT)[1,"Pr(>F)"]
    }
  }
  if((nNoBin==0) & ("AUC" %in% outputs)){
    regressionAUC <- lm(res$AUCv.global[(res$AUCv.global$variable=="Mean" & res$AUCv.global$dataP %in% (dimlabP[2:(npermut+1)])),"mean"] ~ res$prctGlob.Ychange.values[2:(npermut+1),2])
    res$reg.GlobalRes_prctYchange["AUCvGlobal","coeff"]  <- regressionAUC$coefficients[2]
    res$reg.GlobalRes_prctYchange["AUCvGlobal","pvalue"] <- anova(regressionAUC)[1,"Pr(>F)"]
  }
  if(("ER" %in% outputs) | ("AUC" %in% outputs)){
    res$reg.GlobalRes_prctYchange <- data.frame(RegInFunction.prctYchange = rownames(res$reg.GlobalRes_prctYchange),res$reg.GlobalRes_prctYchange)
    rownames(res$reg.GlobalRes_prctYchange) <- NULL
  }
  
  # student test
  
  if(is.null(nbObsPermut)!=TRUE){
    if(("ER" %in% outputs) | ("AUC" %in% outputs)){
      res$ttestMeanERv <- matrix(NA,nrow = 4, ncol = 3, dimnames=list(c("ttestMeanERvMax","ttestMeanERvGravity","ttestMeanERvThreshold","ttestMeanAUCv"),c("mean.noPermut","mean.dataP","p.value")))
    }
    if("ER" %in% outputs){
      if("max" %in% algo){
        ttestMeanERvMax <- t.test(res$ErrorRateVglobal.max[(res$ErrorRateVglobal.max$variable=="global" & res$ErrorRateVglobal.max$dataP %in% (dimlabP[2:(npermut+1)])),"mean"], 
                                  mu=res$ErrorRateVglobal.max[(res$ErrorRateVglobal.max$variable=="global" & res$ErrorRateVglobal.max$dataP == dimlabP[1]),"mean"])
        res$ttestMeanERv["ttestMeanERvMax",] <- unlist(ttestMeanERvMax[c("null.value","estimate","p.value")])
      }
      if("gravity" %in% algo){
        ttestMeanERvGravity <- t.test(res$ErrorRateVglobal.gravity[(res$ErrorRateVglobal.gravity$variable=="global" & res$ErrorRateVglobal.gravity$dataP %in% (dimlabP[2:(npermut+1)])),"mean"], 
                                      mu=res$ErrorRateVglobal.gravity[(res$ErrorRateVglobal.gravity$variable=="global" & res$ErrorRateVglobal.gravity$dataP == dimlabP[1]),"mean"])
        res$ttestMeanERv["ttestMeanERvGravity",] <- unlist(ttestMeanERvGravity[c("null.value","estimate","p.value")])
      }
      if("threshold" %in% algo){
        ttestMeanERvThreshold <- t.test(res$ErrorRateVglobal.threshold[(res$ErrorRateVglobal.threshold$variable=="global" & res$ErrorRateVglobal.threshold$dataP %in% (dimlabP[2:(npermut+1)])),"mean"], 
                                        mu=res$ErrorRateVglobal.threshold[(res$ErrorRateVglobal.threshold$variable=="global" & res$ErrorRateVglobal.threshold$dataP == dimlabP[1]),"mean"])
        res$ttestMeanERv["ttestMeanERvThreshold",] <- unlist(ttestMeanERvThreshold[c("null.value","estimate","p.value")])
      }
    }
    if((nNoBin==0) & ("AUC" %in% outputs)){
      ttestMeanAUCv <- t.test(res$AUCv.global[(res$AUCv.global$variable=="Mean" & res$AUCv.global$dataP %in% (dimlabP[2:(npermut+1)])),"mean"], 
                              mu=res$AUCv.global[(res$AUCv.global$variable=="Mean" & res$AUCv.global$dataP == dimlabP[1]),"mean"])
      res$ttestMeanERv["ttestMeanAUCv",] <- unlist(ttestMeanAUCv[c("null.value","estimate","p.value")])
    }
    if(("ER" %in% outputs) | ("AUC" %in% outputs)){
      #res$ttestMeanERv           <- res$ttestMeanERv[complete.cases(res$ttestMeanERv), ]
      res$ttestMeanERv           <- data.frame(ttest= rownames(res$ttestMeanERv),res$ttestMeanERv, stringsAsFactors = TRUE)
      rownames(res$ttestMeanERv) <- NULL
    }
  }
  res$call <- match.call()
  class(res) <- c("permut_mbplsda")
  return(res)
}

