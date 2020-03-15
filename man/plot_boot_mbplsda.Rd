\name{plot_boot_mbplsda}
\alias{plot_boot_mbplsda}

\title{
Plot the results of the fonction boot_mbplsda in a pdf file
}

\description{
Fonction to draw the results of the fonction boot_mbplsda (2-fold cross-validated parameter values) in a pdf file
}

\usage{
plot_boot_mbplsda(obj, filename = "PlotBootstrapMbplsda", propbestvar = 0.5)
}

\arguments{
  \item{obj}{
object type list containing the results of the fonction boot_mbplsda
}
  \item{filename}{
a string of characters indicating the given pdf filename
}
  \item{propbestvar}{
numeric value between 0 and 1, indicating the pourcentage of variables with the best VIPc values to plot
}
}

\details{
no details are needed
}

\value{
no numeric result
}

\references{
Efron, B., Tibshirani, R.J. (1994). An Introduction to the Bootstrap. Chapman and Hall-CRC Monographs on Statistics and Applied Probability, Norwell, Massachusetts, United States.

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at 12emes Journees Scientifiques RFMF, Clermont-Ferrand, FRA(05-21-2019 - 05-23-2019).

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics, 15(10):134

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2020). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at Chimiometrie 2020, Liege, BEL(01-27-2020 - 01-29-2020).
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=boot_mbplsda]{boot_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
\donttest{
data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical, nutrition = nutrition, omics = omics))
disjonctif <- (disjunctive(status))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
resboot <- boot_mbplsda(modelembplsQ, optdim = ncpopt, nrepet = 30, cpus=1)
plot_boot_mbplsda(resboot,"plotBoot_nf1_30rep", propbestvar=0.20)
}
}

\keyword{multivariate}
