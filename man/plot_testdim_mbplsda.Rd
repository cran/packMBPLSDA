\name{plot_testdim_mbplsda}
\alias{plot_testdim_mbplsda}

\title{
Plot the results of the fonction testdim_mbplsda in a pdf file
}

\description{
Fonction to draw the results of the fonction testdim_mbplsda (cross validated prediction error rates, or aera under ROC curve, in function of the number of components in the model) in a pdf file
}

\usage{
plot_testdim_mbplsda(obj, filename = "PlotTestdimMbplsda")
}

\arguments{
  \item{obj}{
object type list containing the results of the fonction testdim_mbplsda
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
Stone, M. (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the Royal Statistical Society B, 36(2), 111-147.

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at 12emes Journees Scientifiques RFMF, Clermont-Ferrand, FRA(05-21-2019 - 05-23-2019).

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics, 15(10):134

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2020). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at Chimiometrie 2020, Liege, BEL(01-27-2020 - 01-29-2020).
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=testdim_mbplsda]{testdim_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
\donttest{
data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical[,1:10], 
nutrition = nutrition[,1:10], omics = omics[,1:20]))
disjonctif <- (disjunctive(status))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
bloYobs <- 2
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 3)
resdim <- testdim_mbplsda(object=modelembplsQ, nrepet = 30, threshold = 0.5, 
bloY=bloYobs, cpus=1, algo = c("max"), outputs = c("ER"))
plot_testdim_mbplsda(resdim, "plotTDim")
}
}

\keyword{multivariate}

