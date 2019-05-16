\name{plot_cvpred_mbplsda}
\alias{plot_cvpred_mbplsda}

\title{
Plot the results of the fonction cvpred_mbplsda in a pdf file
}

\description{
Fonction to draw the results of the fonction cvpred_mbplsda (2-fold cross-validated predictions) in a pdf file
}

\usage{
plot_cvpred_mbplsda(obj, filename = "PlotCVpredMbplsda")
}

\arguments{
  \item{obj}{
object type list containing the results of the fonction cvpred_mbplsda
}
  \item{filename}{
a string of characters indicating the given pdf filename
}
}

\details{
no details are needed
}

\value{
no numeric result
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)

Stone, M. (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the Royal Statistical Society B, 36(2), 111-147.
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=cvpred_mbplsda]{cvpred_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical[,1:10], 
nutrition = nutrition[,1:10], omics = omics[,1:20]))
disjonctif <- (disjunctive(status))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
bloYobs <- 2
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", 
scannf = FALSE, nf = 2)
CVpred <- cvpred_mbplsda(modelembplsQ, nrepet = 30, threshold = 0.5, bloY=bloYobs, 
optdim=ncpopt, cpus = 1, algo = c("max"))
plot_cvpred_mbplsda(CVpred,"plotCVPred_nf1_30rep")

\donttest{
data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical, 
nutrition = nutrition, omics = omics))
disjonctif <- (disjunctive(status))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
bloYobs <- 2
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", 
scannf = FALSE, nf = 2)
CVpred <- cvpred_mbplsda(modelembplsQ, nrepet = 90, threshold = 0.5, bloY=bloYobs, 
optdim=ncpopt, cpus = 1, algo = c("max"))
plot_cvpred_mbplsda(CVpred,"plotCVPred_nf1_90rep")
}
}

\keyword{multivariate}

