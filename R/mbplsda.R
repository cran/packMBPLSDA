mbplsda <- function(dudiY, ktabX, scale = TRUE, option = c("uniform", "none"), scannf = TRUE, nf = 2) {
  
  ## -------------------------------------------------------------------------------
  ##  		       Some tests
  ##--------------------------------------------------------------------------------
  
  if (!inherits(dudiY, "dudi"))
    stop("object 'dudi' expected")
  if (!inherits(ktabX, "ktab"))
    stop("object 'ktab' expected")
  if (any(row.names(ktabX) != row.names(dudiY$tab)))
    stop("ktabX and dudiY must have the same rows")
  if (!(all.equal(ktabX$lw/sum(ktabX$lw), dudiY$lw/sum(dudiY$lw))))
    stop("ktabX and dudiY must have the same row weights")
  if (nrow(dudiY$tab) < 6)
    stop("Minimum six rows are required")
  if (any(ktabX$blo < 2))
    stop("Minimum two variables per explanatory block are required")
  if (!(is.logical(scale)))
    stop("Non convenient selection for scaling")
  if (!(is.logical(scannf)))
    stop("Non convenient selection for scannf")
  if (nf < 0)
    nf <- 2
  
  option <- match.arg(option)
  
  
  ## -------------------------------------------------------------------------------
  ##			Arguments and data transformation
  ## -------------------------------------------------------------------------------
 
  nr           <- nrow(dudiY$tab)
  meanY        <- colMeans(as.matrix(dudiY$tab))
  sdY          <- apply(as.matrix(dudiY$tab), 2, sd)   # biaised sd
  blo          <- ktabX$blo
  nblo         <- length(ktabX$blo)
  meanX        <- colMeans(cbind.data.frame(lapply(unclass(ktabX)[1 : nblo], scale, center = FALSE, scale = FALSE)))
  names(meanX) <- col.names(ktabX)
  sdX          <- apply(cbind.data.frame(lapply(unclass(ktabX)[1 : nblo], scale, center = FALSE, scale = FALSE)), 2, sd) * sqrt((nr-1)/nr)  # biaised sd
 
  
  ## Pre-processing of the data frames (centering by default + biaised variances)
  Y  <- as.matrix(scalewt(dudiY$tab, center = TRUE, scale = scale))
  Xk <- lapply(unclass(ktabX)[1 : nblo], scalewt, wt = ktabX$lw, center = TRUE, scale = scale)
  
  
  ## Block weighting (/ Unbiaised inertia)
  if (option[1] == "uniform"){
    inertia <- lapply(1:nblo, function(k) inertie(Xk[[k]]))
    Xk      <- lapply(1:nblo, function(k) Xk[[k]] / sqrt(inertie(Xk[[k]])))
  }
  X           <- cbind.data.frame(Xk)
  colnames(X) <- col.names(ktabX)
  
  
  ## Parameters
  ncolX   <- ncol(X)
  ncolY   <- ncol(Y)
  if (scannf == TRUE) {maxdim  <- qr(X)$rank}
  if (scannf == FALSE){maxdim  <- min(nf,qr(X)$rank)}
  
  
  # (non-centered) (eventually) reduced and weighted variable means
  meanX1 <- meanX
  meanY1 <- meanY
  if (scale == TRUE){
    meanX1 <- meanX1 / sdX
    meanY1 <- meanY1 / sdY
  } 
  if (option[1] == "uniform"){
    meanX1 <- meanX1 / sqrt(unlist(lapply(1 : nblo, function(k) rep(inertia[[k]], blo[k]))))
  }
  meanX.transf <- meanX1 
  meanY.transf <- meanY1
  
  
  
  ##-----------------------------------------------------------------------
  ##                         Prepare the outputs
  ##-----------------------------------------------------------------------
  
  
  dimlab   <- paste("Ax", 1:maxdim, sep = "")
  
  res      <- list(tabX = X, tabY = as.data.frame(Y), nf = nf, lw = ktabX$lw, X.cw = ktabX$cw, blo = ktabX$blo, rank = maxdim, eig = rep(0, maxdim), TL = ktabX$TL, TC = ktabX$TC)
  res$Yc1  <- matrix(0, nrow = ncolY, ncol = maxdim, dimnames = list(colnames(dudiY$tab), dimlab))
  res$lX   <- res$lY <- matrix(0, nrow = nr, ncol = maxdim, dimnames = list(row.names(dudiY$tab), dimlab))
  res$cov2 <- Ak <- matrix(0, nrow = nblo, ncol = maxdim, dimnames = list(names(ktabX$blo), dimlab))
  res$Tc1  <- lapply(1:nblo, function(k)  matrix(0, nrow = ncol(Xk[[k]]), ncol = maxdim, dimnames = list(colnames(Xk[[k]]), dimlab)))
  res$TlX  <- rep(list(matrix(0, nrow = nr, ncol = maxdim, dimnames = list(row.names(dudiY$tab), dimlab))), nblo)
  res$faX  <- matrix(0, nrow = ncolX, ncol = maxdim, dimnames = list(col.names(ktabX), dimlab))
  lX1      <- res$lX
  W        <- res$faX
  
  
  
  ##-----------------------------------------------------------------------
  ##     Compute components and loadings by an iterative algorithm
  ##-----------------------------------------------------------------------
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  f1 <- function(x) crossprod(x * res$lw, Y)
  
  for(h in 1 : maxdim) {
    

    ## Compute the matrix M for the eigenanalysis
    M <- lapply(lapply(Xk, f1), crossprod)
    M <- Reduce("+", M)
    
    ## Compute the loadings V and the components U (Y dataset)
    eig.M <- eigen(M)
    
    if (eig.M$values[1] < sqrt(.Machine$double.eps)) {
      res$rank <- h-1 ## update the rank
      break
    }
    
    res$eig[h] <- eig.M$values[1]    
    res$Yc1[, h]  <- eig.M$vectors[, 1, drop = FALSE]
    res$lY[, h]  <- Y %*% res$Yc1[, h]
    
    ## Compute the loadings Wk and the components Tk (Xk datasets)
    covutk <- rep(0, nblo)   
    for (k in 1 : nblo) {
      res$Tc1[[k]][, h] <- crossprod(Xk[[k]] * res$lw, res$lY[, h])
      res$Tc1[[k]][, h] <- res$Tc1[[k]][, h] / sqrt(sum(res$Tc1[[k]][, h]^2))
      res$TlX[[k]][, h] <- Xk[[k]] %*% res$Tc1[[k]][, h]
      
      covutk[k] <- crossprod(res$lY[, h] * res$lw, res$TlX[[k]][, h])
      res$cov2[k, h] <- covutk[k]^2  # equivalent to apply(crossprod(as.matrix(modelembplsQ$tabY* modelembplsQ$lw), (modelembplsQ$TlX)[[k]][, h])^2,sum, MARGIN = 2)
    }
    
    for(k in 1 : nblo) {
      Ak[k, h]     <- covutk[k] / sqrt(sum(res$cov2[,h]))
      res$lX[, h]  <- res$lX[, h] + Ak[k, h] * res$TlX[[k]][, h]
    }
    
    lX1[, h] <- res$lX[, h] / sqrt(sum(res$lX[, h]^2))
    W[, h]   <- tcrossprod(MASS::ginv(crossprod(X)), X) %*% res$lX[, h]    ## NB : use ginv to avoid NA in coefficients (collinear system)
    
    ## Deflation of the Xk datasets on the global components T
    Xk <- lapply(Xk, function(y) lm.wfit(x = as.matrix(res$lX[, h]), y = y, w = res$lw)$residuals)
    X  <- as.matrix(cbind.data.frame(Xk))
  }
  
  
  
  ##-----------------------------------------------------------------------
  ##     Compute regressions coefficients
  ##-----------------------------------------------------------------------
  
  ## Use of the original (and not the deflated) datasets X and Y	
  X <- as.matrix(res$tabX)
  Y <- as.matrix(res$tabY)
  
  
  ## Computing the regression coefficients of X onto the global components T (Wstar)
  res$faX[, 1] <- W[, 1, drop = FALSE]
  A           <- diag(ncolX)
  Xd          <- X
  if(maxdim >= 2){
    for(h in 2 : maxdim){
      a            <- crossprod(lX1[, h-1], Xd) / sqrt(sum(res$lX[, h-1]^2))
      A            <- A %*% (diag(ncolX) - W[, h-1] %*% a)
      res$faX[, h] <- A %*% W[, h]
      Xd           <- Xd - tcrossprod(lX1[, h-1]) %*% Xd
    }
  }
  
  
  ##  Computing the regression coefficients of X onto Y (Beta) for X and Y centred and eventually reduced and weighted
  res$Yco <-  t(Y) %*% diag(res$lw) %*% res$lX
  norm.li <- diag(crossprod(res$lX * sqrt(res$lw)))
  if (maxdim == 1){
    res$XYcoef <- lapply(1:ncolY, function(q) as.matrix(apply(sweep(res$faX, 2 , res$Yco[q, ] / norm.li, "*"), 1, cumsum)))
  } else {
    res$XYcoef <- lapply(1:ncolY, function(q) t(apply(sweep(res$faX, 2 , res$Yco[q, ] / norm.li, "*"), 1, cumsum)))
  }
  names(res$XYcoef) <- colnames(dudiY$tab)
  
  
  ## Correct the regression coefficients in case of block weighting (option == uniform)
  beta1 <- res$XYcoef
  if (option[1] == "uniform"){
    index  <- as.factor(rep(1:nblo, blo))
    for (q in 1 : ncolY){
      if(maxdim==1){
        beta1[[q]] <- cbind(unlist(lapply(1:nblo, function(k) beta1[[q]][index == k, ] / sqrt(inertia[[k]]))))
        #beta1[[q]] <- do.call(rbind, lapply(1:nblo, function(k) beta1[[q]][index == k, ] / sqrt(inertia[[k]])))
      } else {
        beta1[[q]] <- do.call(rbind, lapply(1:nblo, function(k) beta1[[q]][index == k, ] / sqrt(inertia[[k]])))
      }
    }
  }
  
  
  
  ## Correct the regression coefficients in case of scaling (scale = TRUE)
  if (scale == TRUE){
    beta1 <- lapply(1:ncolY, function(q) beta1[[q]] * matrix(rep(sdY[q], each=ncolX * maxdim), ncol = maxdim) / t(matrix(rep(sdX, each=maxdim), nrow=maxdim)))
  }
  res$XYcoef.raw <- beta1
  colnames(res$XYcoef.raw) <- colnames(res$XYcoef)

  
  ## Compute the intercept
  res$intercept        <- lapply(1:ncolY, function(q) (meanY.transf[q] - meanX.transf %*% res$XYcoef[[q]]))
  res$intercept.raw    <- lapply(1:ncolY, function(q) (meanY[q] - meanX %*% res$XYcoef.raw[[q]]))
  names(res$intercept) <- names(res$intercept.raw) <- colnames(dudiY$tab) 
  
     
  
  
  ##-----------------------------------------------------------------------
  ##   		Variable and block importances
  ##-----------------------------------------------------------------------
  
  ## Block importances
  res$bip <- Ak^2
  if (nblo == 1 | res$rank ==1){
    res$bipc <- res$bip
  } else {
    res$bipc <- t(sweep(apply(sweep(res$bip, 2, res$eig, "*") , 1, cumsum), 1, cumsum(res$eig), "/"))
  }
    
  
  ## Variable importances
  WcarreAk <- res$faX^2 * res$bip[rep(1:nblo, ktabX$blo),]
  res$vip  <- sweep(WcarreAk, 2, colSums(WcarreAk), "/")
  if (res$rank ==1){ ## ATTENTION MODIFIED in mbpls ade4: if (nblo == 1 | res$rank ==1)
    res$vipc <- res$vip
  } else {
    res$vipc <- t(sweep(apply(sweep(res$vip, 2, res$eig, "*") , 1, cumsum), 1, cumsum(res$eig), "/"))
    ## ATTENTION in Mixomics and other packages VIP=res$vipc
    # cor2 = as.matrix(cor(Y$tab, res$lX, use="pairwise")^2, nrow=q)
    # for h = 1, VIP = faX[,1]^2
    # for h > 1, VIP = apply(cor2[, 1:h], 2, sum) %*% t(faX[,h]^2)/ sum(apply(cor2[, 1:h], 2, sum))
  }
    
  
  
  
  ##-----------------------------------------------------------------------
  ##			         Modify the outputs
  ##-----------------------------------------------------------------------
  
  if (scannf == TRUE){
    barplot(res$eig[1:res$rank])
    cat("Select the number of global components: ")
    res$nf <- as.integer(readLines(n = 1))
  }
  
  if(res$nf > res$rank)
    res$nf <- res$rank
  
  ## keep results for the nf dimensions (except eigenvalues and lX)
  res$eig <- res$eig[1:res$rank]
  res$lX  <- data.frame(res$lX[, 1:res$rank])
  colnames(res$lX) <- dimlab[1:res$rank]
  res$Tc1 <- do.call("rbind", res$Tc1)
  res$TlX <- do.call("rbind", res$TlX)
  
  res           <- modifyList(res, lapply(res[c("Yc1", "Yco", "lY", "Tc1", "TlX", "cov2", "faX", "vip", "vipc", "bip", "bipc")], function(x) x[, 1:res$nf, drop = FALSE]))
  res$XYcoef    <- lapply(res$XYcoef, function(x) x[, 1:res$nf, drop = FALSE])
  res$intercept <- lapply(res$intercept, function(x) x[, 1:res$nf, drop = FALSE])
  res$call      <- match.call()
  class(res)    <- c("mbplsda")
  return(res)
}



# ---------------------------------------------------------------------------
# 5. Additional functions
# ---------------------------------------------------------------------------


# Generalized inverse of a matrix (from MASS)
ginv <- function(X, tol = sqrt(.Machine$double.eps)){
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
#  Xsvd <- svd(X,LINPACK = TRUE)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}


# Matrix (biaised) inertia
inertie <- function(tab){
  tab <- as.matrix(scale(tab, center = TRUE, scale = FALSE))
  V   <- tab %*% t(tab) / (dim(tab)[1])
  sum(diag(V))
}
