\name{pred_mbplsda}
\alias{pred_mbplsda}

\title{
Observed parameters and predicted categories from a multi-block partial least squares discriminant model
}

\description{
Fonction to perform categories predictions from a multi-block partial least squares discriminant model.
}

\usage{
pred_mbplsda(object, optdim , threshold = 0.5, bloY, 
algo = c("max", "gravity", "threshold"))
}

\arguments{
  \item{object}{
an object created by mbplsda
}
  \item{optdim}{
integer indicating the (optimal) number of components of the multi-block partial least squares discriminant model
}
  \item{threshold}{
numeric indicating the threshold, between 0 and 1, to consider the categories are predicted with the threshold prediction method.
}
  \item{bloY}{
integer vector indicating the number of categories per variable of the Y-block.
}
  \item{algo}{
character vector indicating the method(s) of prediction to use (see details)
}
}

\details{
Three different algorithms are available to predict the categories of observations.
In the max, and respectively the threshold algorithms, numeric values are calculated from the matrix of explanatory variables and the regression coefficients. Then, the predicted categorie for each variable of the Y-block is the one which corresponds to the higher predicted value, respectively to the values higher than the indicated threshold.
In the gravity algorithm, predicted scores of the observations on the components are calculated. Then, each observation is assigned to the observed category of which it is closest to the barycentre in the component space.
}

\value{
	\item{XYcoef}{list of matrices of the regression coefficients of the whole explanatory dataset onto the dependent dataset}
	\item{VIPc}{cumulated variable importances for a given number of dimensions}
	\item{BIPc}{cumulated block importances for a given number of dimensions}
	\item{faX}{matrix containing the global variable loadings associated with the global explanatory dataset}
	\item{lX}{matrix of the global components associated with the whole explanatory dataset(scores of the individuals)}
	\item{ConfMat.ErrorRate}{confidence matrix and prediction error rate per category}
	\item{ErrorRate.global}{confidence matrix and prediction error rate, per Y-block variable and overall}
	\item{PredY.max}{predictions and accuracy of predictions with the "max" algorithm}
	\item{PredY.gravity}{predictions and accuracy of predictions with the "gravity" algorithm}
	\item{PredY.threshold}{predictions and accuracy of predictions with the "threshold" algorithm}
	\item{AUC}{aera under ROC cuve value and 95\% confidence interval, per category, per Y-block variable and overall}
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=plot_pred_mbplsda]{plot_pred_mbplsda}}
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
bloYobs <- 2
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
predictions <- pred_mbplsda(modelembplsQ, optdim = ncpopt, threshold = 0.5, bloY=bloYobs, 
algo = c("max", "gravity", "threshold"))
}

\keyword{multivariate}