\name{plot_permut_mbplsda}
\alias{plot_permut_mbplsda}

\title{
Plot the results of the fonction permut_mbplsda in a pdf file
}

\description{
Fonction to draw the results of the fonction permut_mbplsda (plot and regression line of cross validated prediction error rates, evaluated on the validation datasets, in function of the percent of modified Y-block values) in a pdf file
}

\usage{
plot_permut_mbplsda(obj, filename = "PlotPermutationTest", 
MainPlot = "Permutation test results \n (subset of validation)")
}

\arguments{
  \item{obj}{
object type list containing the results of the fonction permut_mbplsda 
}
  \item{filename}{
a string of characters indicating the given pdf filename
}
  \item{MainPlot}{
a string of characters indicating the given main title
}
}

\details{
no details are needed
}

\value{
no numeric result
}

\references{
Westerhuis, J.A., Hoefsloot, H.C.J., Smit, S., Vis, D.J., Smilde, A.K., van Velzen, E.J.J., van Duijnhoven, J.P.M., van Dorsten, F.A. (2008). Assessment of PLSDA cross validation. Metabolomics, 4, 81-89.

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at 12emes Journees Scientifiques RFMF, Clermont-Ferrand, FRA(05-21-2019 - 05-23-2019).

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics, 15(10):134

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2020). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at Chimiometrie 2020, Liege, BEL(01-27-2020 - 01-29-2020).
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=permut_mbplsda]{permut_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
\donttest{
data(status)
data(medical)
data(omics)
ktabX <- ktab.list.df(list(medical = medical[1:17,1:2], omics = omics[1:17,1:2]))
disjonctif <- (disjunctive(data.frame(status=status[1:17,], 
row.names = rownames(status)[1:17])))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
bloYobs <- 2
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 1)
ncpopt <- 1
rtsPermut <- permut_mbplsda(modelembplsQ, nrepet = 3, npermut = 3, optdim = ncpopt, 
outputs = c("ER"), bloY=bloYobs, nbObsPermut = 10, cpus = 1, algo = c("max"))
plot_permut_mbplsda(rtsPermut,"plotPermut_nf1_3rep_3perm")
}

\donttest{
data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical[1:20,], omics = omics[1:20,]))
disjonctif <- (disjunctive(data.frame(status=status[1:20,], 
row.names = rownames(status)[1:20])))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
bloYobs <- 2
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 1)
ncpopt <- 1
rtsPermut <- permut_mbplsda(modelembplsQ, nrepet = 30, npermut = 100, optdim = ncpopt, 
outputs = c("ER"), bloY=bloYobs, nbObsPermut = 10, cpus = 1, algo = c("max"))
plot_permut_mbplsda(rtsPermut,"plotPermut_nf1_30rep_100perm")
}
}

\keyword{multivariate}

