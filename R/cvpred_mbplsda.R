
# -----------------------------------------------------------------------------------------
# cvpred.multiblock plsda
# -----------------------------------------------------------------------------------------


cvpred_mbplsda <- function(object, nrepet = 100, threshold = 0.5, bloY, optdim, cpus=1, algo=c("max","gravity","threshold")){

# bloY = NEEDED VECTOR = nb of categories by Y block variable
  
  if (!inherits(object, "mbplsda")) 
    stop("Object of type 'mbplsda' expected")
  
  # step 0. packages 
  
#  library(parallel)  # jobs repartition 
#  library(doParallel) # iterations 
#  library(foreach)
#  library(pROC) # auc

 
  
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
  Mod    <- unlist(sapply(1:Ky, function(x) rep(c(1:bloY[x]))))
  cnames <- colnames(Y$tab) #paste0("Var",Var,"Mod",Mod)
  
  nblo   <- length(X$blo)  # nb X blocks
  blo    <- sapply(1:nblo, function(k) dim(X[[k]])[2]) # nb variables by X block
  Bloc   <- as.factor(rep(1:nblo, blo))
  
  nNoBin <- sum(bloY!=2)  # nb no binary Y variables
  
  dimlab <- paste("Ax", optdim, sep = "")
  
  # step 2. Preparation of the parallelized processing
  #nodes <- detectCores()
  cl    <- makeCluster(cpus, type="PSOCK") # initialisation 
  registerDoParallel(cl) # cl scripts creation
  on.exit(stopCluster(cl))

  # step 3. strat of parallelisation
  resForeach <- foreach(i = 1:nrepet, .export=c("mbplsda", "inertie", "ginv"), .packages=c("ade4","pROC"), .errorhandling="remove") %dopar%{

    set.seed(seed=i)
    
    ## I.1. Dividing X and Y into calibration (Xc, Yc) and validation (Xv, Yv) datasets. 
    s  <- sample(x = nr, size = Nc)  
    Xc <- X[, s, ]
    Xv <- X[, -s, ]
    Yc <- Y[s, ]
    Yv <- Y[-s, ]
    
    rnamesXc <- row.names(Xc)
    rnamesXv <- row.names(Xv)
    rnamesYc <- row.names(Yc$tab)
    rnamesYv <- row.names(Yv$tab)
    
    ## nb "1" by Y category
    nbY1c <- sapply(1:q, function(g) sum(Yc$tab[, g] == 1))
    nbY1v <- sapply(1:q, function(g) sum(Yv$tab[, g] == 1))
    
    
    ## I.2. Application of mbplsda on calibration/validation datasets       
    rescal <- do.call(method, list(dudiY = Yc, ktabX = Xc, scale = scale, option = option, scannf = FALSE, nf = optdim))
    resval <- do.call(method, list(dudiY = Yv, ktabX = Xv, scale = scale, option = option, scannf = FALSE, nf = optdim))                   
    
    ## maximum nb of components
    #H = min(rescal$rank, resval$rank, h)
    
    
    ## I.3. variable matrix
    ## raw Xc  
    Xc.mat <- cbind.data.frame(unclass(Xc)[1:nblo], stringsAsFactors = TRUE)
    
    ## means and biased sd of Xc variables (to use with mbpls ade4)
    rescal$meanX        <- colMeans(Xc.mat) # equivalent to rescal$meanX
    rescal$sdX          <- apply(Xc.mat, 2, sd) * sqrt((Nc-1)/Nc)  # biaised sd, equivalent to rescal$sdX
    
    ## Xv: raw Xv , Xv centred with Xc means, Xv centred reduced weighted with means, sd, inertia of Xc blocks
    Xv.raw <- cbind.data.frame(unclass(Xv)[1:nblo], stringsAsFactors = TRUE)
    
    Xv.c   <- sweep(Xv.raw, 2, rescal$meanX, FUN="-") # equivalent to Xv.raw - rep(rescal$meanX, each=Nv)
    
    if(scale==TRUE){
      Xv.cr  <- sweep(Xv.c, 2, rescal$sdX, FUN="/") # equivalent to Xv.c / rep(rescal$sdX, each=Nv)
      if(option=="uniform"){
        Xv.crw <- Xv.cr * sqrt(matrix(rep(rescal$X.cw, each=Nv), nrow=Nv)) # rescal$X.cw = column weighs
      }else{
        Xv.crw <- Xv.cr
      }
    }
    
    if(scale==FALSE){
      Xv.cr  <- Xv.c
      if(option=="uniform"){
        Xv.crw <- Xv.cr * sqrt(matrix(rep(rescal$X.cw, each=Nv), nrow=Nv)) # rescal$X.cw = column weighs
      }else{
        Xv.crw <- Xv.cr
      }
    }

    # J.1. outputs preparation for optdim

    ## index of predicted categories
    classYc           <- matrix(NA,Nc,ncol=3*Ky,dimnames=list(rnamesYc, c(paste0("Ymax_Var",1:Ky),paste0("Ycentroid_Var",1:Ky), paste0("Ythreshold_Var",1:Ky))))
    classYv           <- matrix(NA,Nv,ncol=3*Ky,dimnames=list(rnamesYv, c(paste0("Ymax_Var",1:Ky),paste0("Ycentroid_Var",1:Ky), paste0("Ythreshold_Var",1:Ky))))

    if("max" %in% algo){
      ## predicted disjunctive table with max
      ClassPredYc.max <- matrix(0,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
      ClassPredYv.max <- matrix(0,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
    }
    
    if("gravity" %in% algo){
      ## predicted disjunctive table with gravity
      ClassPredYc.gravity <- matrix(0,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
      ClassPredYv.gravity <- matrix(0,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
    }
    
    if("threshold" %in% algo){
      ## predicted disjunctive table with threshold 
      ClassPredYc.threshold <- matrix(NA,Nc,ncol=q,dimnames=list(rnamesYc,cnames))
      ClassPredYv.threshold <- matrix(NA,Nv,ncol=q,dimnames=list(rnamesYv,cnames))
    }
    
    
    ## error rates by variable and overall
#      ErrorRateGlobal    <- matrix(NA,nrow=6,ncol=(Ky+1),
#                                   dimnames=list(c("ErrRateCmaxGlobal","ErrRateVmaxGlobal","ErrRateCgravityGlobal","ErrRateVgravityGlobal","ErrRateCthresholdGlobal","ErrRateVthresholdGlobal"), c(paste0("Var",1:Ky),"global")))
    
    
    # J.2. coefficients for raw data
    XYcoef.raw.cal    <- sapply(rescal$XYcoef.raw, function(x) x[, optdim])
    intercept.raw.cal <- sapply(rescal$intercept.raw, function(x) x[, optdim])
    

    # J.3. Yc and Yv predictions by optdim
    predYc        <- matrix(rep(intercept.raw.cal, each=Nc), ncol=q) + as.matrix(Xc.mat) %*% XYcoef.raw.cal
    predYv        <- matrix(rep(intercept.raw.cal, each=Nv), ncol=q) + as.matrix(Xv.raw) %*% XYcoef.raw.cal
    colnames(predYc) <- colnames(predYv) <- cnames

    
    # J.4. RMSE
    residYc       <- as.matrix(Yc$tab) - predYc     
    RMSEC         <- sqrt(sum(residYc^2) / (Nc * q))
    residYv       <- as.matrix(Yv$tab) - predYv
    RMSEV         <- sqrt(sum(residYv^2) / (Nv * q))
    
    
    # J.5. predictions and error rates by category with max
    if("max" %in% algo){
      for (k in 1:Ky){
        
        ## index of predicted categories
        classYc[,k] <- sapply(1:Nc, function(n) which.max(predYc[n,Var == k]))
        classYv[,k] <- sapply(1:Nv, function(n) which.max(predYv[n,Var == k]))
        
        ## predicted disjunctive table
        for(n in 1:Nc){
          if(k==1) (ClassPredYc.max[n,classYc[n,1]] <- 1)
          if(k>1)  (ClassPredYc.max[n,(sum(bloY[1:(k-1)])+classYc[n,k])] <- 1)
        }
        for(n in 1:Nv){
          if(k==1) (ClassPredYv.max[n,classYv[n,1]] <- 1)
          if(k>1)  (ClassPredYv.max[n,(sum(bloY[1:(k-1)])+classYv[n,k])] <- 1)
        }
      }
    }
    

    
    # J.6. predictions and error rates by category with gravity
    if("gravity" %in% algo){
      ## group barycenters on global components
      
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
      
      ## dist to barycenters
      dist.eucl.gravity.Yc    <- sapply(1:q, function(g) apply((rescal$lX[,1:optdim] -  matrix(rep(Gravity[g, ], each = Nc), nrow = Nc))^2, 1, sum))
      dist.eucl.gravity.Yv    <- sapply(1:q, function(g) apply(((as.matrix(Xv.crw)%*%rescal$faX[,1:optdim]) -  matrix(rep(Gravity[g, ], each = Nv), nrow = Nv))^2, 1, sum))
      
      for (k in 1:Ky){
        
        ## index of predicted categories
        classYc[,(k+Ky)] <- sapply(1:Nc, function(n) which.min(dist.eucl.gravity.Yc[n,Var == k]))
        classYv[,(k+Ky)] <- sapply(1:Nv, function(n) which.min(dist.eucl.gravity.Yv[n,Var == k]))
        
        ## predicted disjunctive table
        for(n in 1:Nc){
          if(k==1) (ClassPredYc.gravity[n,classYc[n,(1+Ky)]] <- 1)
          if(k>1)  (ClassPredYc.gravity[n,(sum(bloY[1:(k-1)])+classYc[n,(k+Ky)])] <- 1)
        }
        for(n in 1:Nv){
          if(k==1) (ClassPredYv.gravity[n,classYv[n,(1+Ky)]] <- 1)
          if(k>1)  (ClassPredYv.gravity[n,(sum(bloY[1:(k-1)])+classYv[n,(k+Ky)])] <- 1)
        }
      }
      
    }
    
  
    
    # J.7. predictions and error rates by category with threshold
    if("threshold" %in% algo){
      ## predicted disjunctive table
      ClassPredYc.threshold[predYc>=threshold] <- 1
      ClassPredYv.threshold[predYv>=threshold] <- 1
      ClassPredYc.threshold[predYc<threshold]  <- 0
      ClassPredYv.threshold[predYv<threshold]  <- 0
      
      for (k in 1:Ky){
        ## index of predicted categories
        classYc[,(2*Ky+k)] <- sapply(1:Nc, function(n) if(sum(ClassPredYc.threshold[n,Var == k], na.rm=T)==1) (which.max(ClassPredYc.threshold[n,Var == k])))
        classYv[,(2*Ky+k)] <- sapply(1:Nv, function(n) if(sum(ClassPredYv.threshold[n,Var == k], na.rm=T)==1) (which.max(ClassPredYv.threshold[n,Var == k])))      
      }
    }
    
      

  
    ## additional matrix in order to dim of predicted indicators matrix are nr * q
    ## foreach results list
    PredForeach <- list(NULL)
    
    if("max" %in% algo){
      ClassPredYc.max        <- rbind(ClassPredYc.max, matrix(NA, nrow=Nv, ncol=q, dimnames=list(rnamesYv,cnames)))
      PredForeach$PREDcM     <- ClassPredYc.max[sort(rownames(ClassPredYc.max)),]
      ClassPredYv.max        <- rbind(ClassPredYv.max, matrix(NA, nrow=Nc, ncol=q, dimnames=list(rnamesYc,cnames)))
      PredForeach$PREDvM     <- ClassPredYv.max[sort(rownames(ClassPredYv.max)),]
    }
     
    if("gravity" %in% algo){ 
      ClassPredYc.gravity    <- rbind(ClassPredYc.gravity, matrix(NA, nrow=Nv, ncol=q, dimnames=list(rnamesYv,cnames)))
      PredForeach$PREDcG     <- ClassPredYc.gravity[sort(rownames(ClassPredYc.gravity)),]
      ClassPredYv.gravity    <- rbind(ClassPredYv.gravity, matrix(NA, nrow=Nc, ncol=q, dimnames=list(rnamesYc,cnames)))
      PredForeach$PREDvG     <- ClassPredYv.gravity[sort(rownames(ClassPredYv.gravity)),]
    }
      
    if("threshold" %in% algo){
      ClassPredYc.threshold  <- rbind(ClassPredYc.threshold, matrix(NA, nrow=Nv, ncol=q, dimnames=list(rnamesYv,cnames)))
      PredForeach$PREDcT     <- ClassPredYc.threshold[sort(rownames(ClassPredYc.threshold)),]
      ClassPredYv.threshold  <- rbind(ClassPredYv.threshold, matrix(NA, nrow=Nc, ncol=q, dimnames=list(rnamesYc,cnames)))
      PredForeach$PREDvT     <- ClassPredYv.threshold[sort(rownames(ClassPredYv.threshold)),]
      
    }
    
    ## I.5. foreach results list
    PredForeach 
  }  
  
  # step 4. End of the parallelized loops

  stopCluster(cl)
  on.exit(stopCluster)
#  resForeach
  
  res <- NULL
  nrepetFE       <- length(resForeach)
  if((nrepetFE<1.5)|(is.null(nrepetFE)==TRUE)){
    stop("No adjustement of models")
  }
  res$TRUEnrepet <- nrepetFE
  
  # step 5. prepare outputs 
 
   ## concatenated results
  if("max" %in% algo){
    PREDcMall <- PREDvMall <- matrix(NA, nrow=nr, ncol=q*nrepetFE, dimnames=list(row.names(Y$tab),paste0(cnames,"_rep",rep(1:nrepetFE,rep(q,nrepetFE)))))
    PREDcM <- PREDvM <- list()
    matPREDcM <- matPREDvM <- matrix(NA, nrow = nr, ncol = q, dimnames = list(row.names(Y$tab), c(paste0("Pred_",cnames))))
    res$matPredYc.max <- res$matPredYv.max <- matrix(NA,nrow=nr, ncol=(q+Ky+1), dimnames=list(rownames(Y$tab), c(paste0("Accuracy_",cnames),paste0("Accuracy_Var",1:Ky),"GlobalAccuracy")))
  }
  if("gravity" %in% algo){
    PREDcGall <- PREDvGall <- matrix(NA, nrow=nr, ncol=q*nrepetFE, dimnames=list(row.names(Y$tab),paste0(cnames,"_rep",rep(1:nrepetFE,rep(q,nrepetFE)))))
    PREDcG <- PREDvG <- list()
    matPREDcG <- matPREDvG <- matrix(NA, nrow = nr, ncol = q, dimnames = list(row.names(Y$tab), c(paste0("Pred_",cnames))))
    res$matPredYc.gravity <- res$matPredYv.gravity <- matrix(NA,nrow=nr, ncol=(q+Ky+1), dimnames=list(rownames(Y$tab), c(paste0("Accuracy_",cnames),paste0("Accuracy_Var",1:Ky),"GlobalAccuracy")))
  }
  if("threshold" %in% algo){
    PREDcTall <- PREDvTall <- matrix(NA, nrow=nr, ncol=q*nrepetFE, dimnames=list(row.names(Y$tab),paste0(cnames,"_rep",rep(1:nrepetFE,rep(q,nrepetFE)))))
    PREDcT <- PREDvT <- list()
    matPREDcT <- matPREDvT <- matrix(NA, nrow = nr, ncol = q, dimnames = list(row.names(Y$tab), c(paste0("Pred_",cnames))))
    res$matPredYc.threshold <- res$matPredYv.threshold <- matrix(NA,nrow=nr, ncol=(q+Ky+1), dimnames=list(rownames(Y$tab), c(paste0("Accuracy_",cnames),paste0("Accuracy_Var",1:Ky),"GlobalAccuracy")))
  }

  
  # step 6. means onrepetitions
  ### concatenated results
  repetitions <- as.factor(rep(1:nrepetFE,rep(q,nrepetFE)))
  levels(repetitions)
  
  
  for(i in 1:nrepetFE){
    if("max" %in% algo){
      PREDcMall[1:(dim(resForeach[[i]][["PREDcM"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDcM"]])
      PREDvMall[1:(dim(resForeach[[i]][["PREDvM"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDvM"]])
    }
    if("gravity" %in% algo){
      PREDcGall[1:(dim(resForeach[[i]][["PREDcG"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDcG"]])
      PREDvGall[1:(dim(resForeach[[i]][["PREDvG"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDvG"]])
    }
    if("threshold" %in% algo){
      PREDcTall[1:(dim(resForeach[[i]][["PREDcT"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDcT"]])
      PREDvTall[1:(dim(resForeach[[i]][["PREDvT"]])[1]),repetitions==i] <- (resForeach[[i]][["PREDvT"]])
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
  
  
  stat.desc   <- function (x){ # x is a matrix, columns are variables, lines are repetitions
    nombre    <- colSums(!is.na(x))
    mod       <- round(colMeans(x,na.rm = TRUE))
    moy       <- round(colMeans(x,na.rm = TRUE),5)
    etype     <- round(apply(x,2,sd, na.rm=TRUE),5)
    quartiles <- round(t(apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),5)
    IC        <- t(apply(x, 2, IC95))
    result    <- cbind.data.frame(nombre, mod, moy, etype, IC, quartiles, stringsAsFactors = TRUE)
    colnames(result) <- c("nb", "ModalValue", "Proba.be1", "sd", "95CIinf", "95CIsup","Q2.5", "median", "Q97.5")
    rownames(result) <- colnames(x)
    return(result)
  }

  ### means and sd results on the nrepetFE repetitions

  if("max" %in% algo){
    for (l in 1:q){
      index <- seq(from=l, to=q*nrepetFE, by=q)
      PREDcM[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDcMall[,index])))
      PREDvM[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDvMall[,index])))
      matPREDcM[,l]  <- round(colMeans(t(PREDcMall[,index]),na.rm = TRUE),0)
      matPREDvM[,l]  <- round(colMeans(t(PREDvMall[,index]),na.rm = TRUE),0)
    }
      
    res$statPredYc.max <- do.call("rbind", PREDcM)
    res$statPredYv.max <- do.call("rbind", PREDvM)
    rownames(res$statPredYc.max) <- rownames(res$statPredYv.max) <- NULL
    
    ## accuracy indicators
    for(n in 1:nr){
      res$matPredYc.max[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDcM[n,Q]-Y$tab[n,Q])^2)))
      res$matPredYv.max[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDvM[n,Q]-Y$tab[n,Q])^2)))
    }
    #res$matPredYv.max[1:nr , 1:q] <- 1-((matPREDvM-Y$tab)^2)
    
    for(n in 1:nr){
      res$matPredYc.max[n,(q+1):(q+Ky)]   <- sapply(1:Ky, function(k)(min(1-(matPREDcM[n,Var == k] - Y$tab[n,Var == k])^2)))
      res$matPredYv.max[n,(q+1):(q+Ky)]   <- sapply(1:Ky, function(k)(min(1-(matPREDvM[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    res$matPredYc.max[,"GlobalAccuracy"]     <- apply(res$matPredYc.max[,1:(q+Ky)], min, MARGIN = 1)
    res$matPredYv.max[,"GlobalAccuracy"]     <- apply(res$matPredYv.max[,1:(q+Ky)], min, MARGIN = 1)
    
    res$matPredYc.max <- data.frame(obs=rownames(Y$tab),Y$tab, matPREDcM, res$matPredYc.max)
    res$matPredYv.max <- data.frame(obs=rownames(Y$tab),Y$tab, matPREDvM, res$matPredYv.max)
    colnames(res$matPredYc.max)[2:(q+1)] <- colnames(res$matPredYv.max)[2:(q+1)] <- paste0("Truth_",cnames)
    rownames(res$matPredYc.max) <- rownames(res$matPredYv.max) <- NULL
  }
  
  if("gravity" %in% algo){
    for (l in 1:q){
      index <- seq(from=l, to=q*nrepetFE, by=q)
      PREDcG[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDcGall[,index])))
      PREDvG[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDvGall[,index])))
      matPREDcG[,l]  <- round(colMeans(t(PREDcGall[,index]),na.rm = TRUE),0)
      matPREDvG[,l]  <- round(colMeans(t(PREDvGall[,index]),na.rm = TRUE),0)
    }
      
    res$statPredYc.gravity <- do.call("rbind", PREDcG)
    res$statPredYv.gravity <- do.call("rbind", PREDvG)
    rownames(res$statPredYc.gravity) <- rownames(res$statPredYv.gravity) <- NULL
    
    ## accuracy indicators
    for(n in 1:nr){
      res$matPredYc.gravity[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDcG[n,Q]-Y$tab[n,Q])^2)))
      res$matPredYv.gravity[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDvG[n,Q]-Y$tab[n,Q])^2)))
    }
    #res$matPredYc.gravity[1:nr , 1:q] <- 1-((matPREDcG-Y$tab)^2)
    #res$matPredYv.gravity[1:nr , 1:q] <- 1-((matPREDvG-Y$tab)^2)
    
    for(n in 1:nr){
      res$matPredYc.gravity[n,(q+1):(q+Ky)]   <- sapply(1:Ky, function(k)(min(1-(matPREDcG[n,Var == k] - Y$tab[n,Var == k])^2)))
      res$matPredYv.gravity[n,(q+1):(q+Ky)]   <- sapply(1:Ky, function(k)(min(1-(matPREDvG[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    res$matPredYc.gravity[,"GlobalAccuracy"]     <- apply(res$matPredYc.gravity[,1:(q+Ky)], min, MARGIN = 1)
    res$matPredYv.gravity[,"GlobalAccuracy"]     <- apply(res$matPredYv.gravity[,1:(q+Ky)], min, MARGIN = 1)
    
    res$matPredYc.gravity <- data.frame(obs=rownames(Y$tab),Y$tab, matPREDcG, res$matPredYc.gravity, stringsAsFactors = TRUE)
    res$matPredYv.gravity <- data.frame(obs=rownames(Y$tab),Y$tab, matPREDvG, res$matPredYv.gravity, stringsAsFactors = TRUE)
    colnames(res$matPredYc.gravity)[2:(q+1)] <- colnames(res$matPredYv.gravity)[2:(q+1)] <- paste0("Truth_",cnames)
    rownames(res$matPredYc.gravity) <- rownames(res$matPredYv.gravity) <- NULL
  }

  if("threshold" %in% algo){
    for (l in 1:q){
      index <- seq(from=l, to=q*nrepetFE, by=q)
      PREDcT[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDcTall[,index])))
      PREDvT[[l]] <- cbind(obs=row.names(Y$tab),cname=rep(cnames[l],nr),stat.desc(t(PREDvTall[,index])))
      matPREDcT[,l]  <- round(colMeans(t(PREDcTall[,index]),na.rm = TRUE),0)
      matPREDvT[,l]  <- round(colMeans(t(PREDvTall[,index]),na.rm = TRUE),0)
    }
      
    res$statPredYc.threshold <- do.call("rbind", PREDcT)
    res$statPredYv.threshold <- do.call("rbind", PREDvT)
    rownames(res$statPredYc.threshold) <- rownames(res$statPredYv.threshold) <- NULL
    
    ## accuracy indicators
    for(n in 1:nr){
      res$matPredYc.threshold[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDcT[n,Q]-Y$tab[n,Q])^2)))
      res$matPredYv.threshold[n , 1:q] <- sapply(1:q, function(Q)(1-((matPREDvT[n,Q]-Y$tab[n,Q])^2)))
    }
    
    for(n in 1:nr){
      res$matPredYc.threshold[n,(q+1):(q+Ky)] <- sapply(1:Ky, function(k)(min(1-(matPREDcT[n,Var == k] - Y$tab[n,Var == k])^2)))
      res$matPredYv.threshold[n,(q+1):(q+Ky)] <- sapply(1:Ky, function(k)(min(1-(matPREDvT[n,Var == k] - Y$tab[n,Var == k])^2)))
    }
    res$matPredYc.threshold[,"GlobalAccuracy"] <- apply(res$matPredYc.threshold[,1:(q+Ky)], min, MARGIN = 1)
    res$matPredYv.threshold[,"GlobalAccuracy"] <- apply(res$matPredYv.threshold[,1:(q+Ky)], min, MARGIN = 1)
    
    res$matPredYc.threshold <- data.frame(obs=rownames(Y$tab), Y$tab, matPREDcT, res$matPredYc.threshold, stringsAsFactors = TRUE)
    res$matPredYv.threshold <- data.frame(obs=rownames(Y$tab), Y$tab, matPREDvT, res$matPredYv.threshold, stringsAsFactors = TRUE)
    colnames(res$matPredYc.threshold)[2:(q+1)] <- colnames(res$matPredYv.threshold)[2:(q+1)] <- paste0("Truth_",cnames)
    rownames(res$matPredYc.threshold) <- rownames(res$matPredYv.threshold) <- NULL
  }

  res$call   <- match.call()
  class(res) <- c("cvpred_mbplsda")
  return(res)
}