\name{boot_mbplsda}
\alias{boot_mbplsda}

\title{
bootstraped simulations for multi-block partial least squares discriminant analysis
}

\description{
Function to perform bootstraped simulations for multi-block partial least squares discriminant analysis, in order to get confidence intervals for regression coefficients, variable loadings, variable and block importances.
}

\usage{
boot_mbplsda(object, nrepet = 199, optdim, cpus = 1, ...)
}

\arguments{
  \item{object}{
an object created by mbplsda
}
  \item{nrepet}{
integer indicating the number of repetitions
}
  \item{optdim}{
integer indicating the optimal number of global components to be introduced in the model
}
  \item{cpus}{
integer indicating the number of cpus to use when running the code in parallel
}
  \item{\dots}{
other arguments to be passed to methods
}
}

\details{
no details are needed
}

\value{
	\item{XYcoef}{mean, standard deviation, quantiles (0.025;0.975), 95\% confidence interval, median for regression coefficients}
	\item{faX}{mean, standard deviation, quantiles (0.025;0.975), 95\% confidence interval, median for variable loadings}
	\item{vipc}{mean, standard deviation, quantiles (0.025;0.975), 95\% confidence interval, median for cumulated variable importances}
	\item{bipc}{mean, standard deviation, quantiles (0.025;0.975), 95\% confidence interval, median for cumulated block importances}
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

\note{
at least 30 bootstrap repetitions may be recommended, more than 100 beeing preferable
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=plot_boot_mbplsda]{plot_boot_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
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
}

\keyword{multivariate}

