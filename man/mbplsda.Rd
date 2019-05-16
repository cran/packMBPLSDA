\name{mbplsda}
\alias{mbplsda}

\title{
Multi-block partial least squares discriminant analysis
}

\description{
Function to perform a multi-block partial least squares discriminant analysis (MBPLSDA) of several explanatory blocks defined as an object of class ktab, to explain a dependent dataset (Y-block) defined as an object of class dudi, in order to get model parameters for the indicated number of components.
}

\usage{
mbplsda(dudiY, ktabX, scale = TRUE, option = c("uniform", "none"), 
scannf = TRUE, nf = 2)
}

\arguments{
  \item{dudiY}{
an object of class dudi containing the dependent variables
}
  \item{ktabX}{
an object of class ktab containing the blocks of explanatory variables
}
  \item{scale}{
logical value indicating whether the explanatory variables should be standardized
}
  \item{option}{
option for the block weighting. If uniform, the weight of each explanatory block is equal to 1/number of explanatory blocks, and the weight of the Y-block is eqyual to 1. If none, the block weight is equal to the block inertia.
}
  \item{scannf}{
logical value indicating whether the eigenvalues bar plot should be displayed
}
  \item{nf}{
integer indicating the number of components to be calculated
}
}

\details{
no details are needed
}

\value{
\item{call}{the matching call}
\item{tabX}{data frame of explanatory variables centered, eventually scaled (if scale=TRUE)and weighted (if option="uniform")}
\item{tabY}{data frame of dependent variables centered, eventually scaled (if scale=TRUE)and weighted (if option="uniform")}
\item{nf}{integer indicating the number of kept dimensions}
\item{lw}{numeric vector of row weights}
\item{X.cw}{numeric vector of column weights for the explanalatory dataset}
\item{blo}{vector of the numbers of variables in each explanatory dataset}
\item{rank}{rank of the analysis}
\item{eig}{numeric vector containing the eigenvalues}
\item{TL}{dataframe useful to manage graphical outputs}
\item{TC}{dataframe useful to manage graphical outputs}
\item{faX}{matrix containing the global variable loadings associated with the global explanatory dataset}
\item{Tc1}{matrix containing the partial variable loadings associated with each explanatory dataset(unit norm)}
\item{Yc1}{matrix of the variable loadings associated with the dependent dataset}
\item{lX}{matrix of the global components associated with the whole explanatory dataset(scores of the individuals)}
\item{TlX}{matrix containing the partial components associated with each explanatory dataset}
\item{lY}{matrix of the components associated with the dependent dataset}
\item{cov2}{squared covariance between lY and TlX}
\item{XYcoef}{list of matrices of the regression coefficients of the whole explanatory dataset onto the dependent dataset}
\item{intercept}{intercept of the regression of the whole explanatory dataset onto the dependent dataset}
\item{XYcoef.raw}{list of matrices of the regression coefficients of the whole raw explanatory dataset onto the raw dependent dataset}
\item{intercept.raw}{intercept of the regression of the whole raw explanatory dataset onto the raw dependent dataset}
\item{bip}{block importances for a given dimension}
\item{bipc}{cumulated block importances for a given number of dimensions}
\item{vip}{variable importances for a given dimension}
\item{vipc}{cumulated variable importances for a given number of dimensions}
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)

Bougeard, S. and Dray, S. (2018) Supervised Multiblock Analysis in R with the ade4 Package.Journal of Statistical Software,86(1), 1-17.
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\note{
This function is coming from the mbpls function of the R package ade4 (application in order to explain a disjunctive table, limitation of the number of calculated components)
}

\seealso{
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
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
}

\keyword{multivariate}

