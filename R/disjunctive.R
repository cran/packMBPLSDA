# -----------------------------------------------------------------------------------------
# needed function to obtain the disjunctive table of the Y block
# -----------------------------------------------------------------------------------------



disjunctive <- function(y){
  binarize <- function(variables){#library(DiscriMiner)
    # binary super-indicator matrix (aka Complete Disjunctive Table)
    # variables: matrix or data.frame with explanatory variables
    
    # make sure variables is a data frame with factors
    fac_check = sapply(variables, class)
    if (!is.data.frame(variables) && any(fac_check != "factor"))
      stop("\n'variables' must be a data frame with factors")
    # no missing values allowed
    if (length(complete.cases(variables)) != nrow(variables))
      stop("\nSorry, no missing values allowed in 'variables'")    
    
    # build super-indicator matrix Z
    my_tdc <- function(X){  
      # Tableau Disjonctive Complete
      # (aka Complete Disjunctive Table)
      # X: data frame with categorical variables as factors
      
      # how many observations
      nobs = nrow(X)
      # how many variables
      nvars = ncol(X)
      # number of categories per variable
      cats_per_var = sapply(X, nlevels)
      # total number of categories
      ncats = sum(cats_per_var)
      
      # build super-indicator matrix Z
      Z = matrix(0, nobs, ncats)
      ini = cumsum(cats_per_var) - cats_per_var + 1
      fin = cumsum(cats_per_var)
      for (j in 1:nvars)
      {
        aux_lev = levels(X[,j])
        aux_mat = matrix(0, nobs, cats_per_var[j])
        for (k in 1:cats_per_var[j])
        {
          tmp <- X[,j] == aux_lev[k]
          aux_mat[tmp,k] = 1
        }
        Z[,ini[j]:fin[j]] = aux_mat
      }
      colnames(Z) = unlist(lapply(X, levels))
      Z
    }
    Z = my_tdc(variables)  
    Z
  }
  
  z <- NULL
  for(i in 1:dim(y)[2]){
    yy <- as.factor(y[,i])
    lev <- levels(yy)
    nlevels <- length(lev)
    zz <- data.frame(y = yy, stringsAsFactors = TRUE)
    z[[i]] <- binarize(zz)
    colnames(z[[i]]) <- paste0(rep(colnames(y)[i],nlevels),"_",colnames(z[[i]]))
  }
  ydisj <- data.frame(lapply (1:dim(y)[2], function(j) (z[[j]])))

  rownames(ydisj) <-rownames(y)
  ydisj
}