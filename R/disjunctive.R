# -----------------------------------------------------------------------------------------
# needed function to obtain the disjunctive table of the Y block
# -----------------------------------------------------------------------------------------
#library(DiscriMiner)


disjunctive <- function(y){
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