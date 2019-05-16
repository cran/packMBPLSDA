\name{cvpred_mbplsda}
\alias{cvpred_mbplsda}

\title{
Cross-validated predicted categories from a multi-block partial least squares discriminant model
}

\description{
Function to perform 2-fold cross-validation for multi-block partial least squares discriminant analysis, in order to get for each observation the cross-validated predicted categories, and the statistical description of the predictions (mean, sd, 95% confidence interval, quantiles (0.025 and 0.975), median) 
}

\usage{
cvpred_mbplsda(object, nrepet = 100, threshold = 0.5, bloY, optdim, cpus = 1, 
algo = c("max", "gravity", "threshold"))
}

\arguments{
  \item{object}{
an object created by mbplsda
}
  \item{nrepet}{
integer indicating the number of repetitions
}
  \item{threshold}{
numeric indicating the threshold, between 0 and 1, to consider the categories are predicted with the threshold prediction method.
}
  \item{bloY}{
integer vector indicating the number of categories per variable of the Y-block.
}
  \item{optdim}{
integer indicating the (optimal) number of components of the multi-block partial least squares discriminant model
}
  \item{cpus}{
integer indicating the number of cpus to use when running the code in parallel
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
\item{TRUEnrepet}{number of repetitions}
\item{matPredYc.max}{with the max algorithm, boolean matrix indicating the cross-validated predicted categories on the calibration datasets, the prediction accuracy for each categorie, each Y-block variable, and overall}
\item{matPredYv.max}{with the max algorithm, boolean matrix indicating the cross-validated predicted categories on the validation datasets, the prediction accuracy for each  categorie, each Y-block variable, and overall}
\item{matPredYc.gravity}{with the gravity algorithm, boolean matrix indicating the cross-validated predicted categories on the calibration datasets, the prediction accuracy for each  categorie, each Y-block variable, and overall}
\item{matPredYv.gravity}{with the gravity algorithm, boolean matrix indicating the cross-validated predicted categories on the validation datasets, the prediction accuracy for each  categorie, each Y-block variable, and overall}
\item{matPredYc.threshold}{with the threshold algorithm, boolean matrix indicating the cross-validated predicted categories on the calibration datasets, the prediction accuracy for each  categorie, each Y-block variable, and overall}
\item{matPredYv.threshold}{with the threshold algorithm, boolean matrix indicating the cross-validated predicted categories on the validation datasets, the prediction accuracy for each  categorie, each Y-block variable, and overall}
\item{statPredYc.max}{with the max algorithm, matrix indicating the statistical description of prediction categories for each observation on the calibration datasets: number of predictions as an observation of the calibration dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
\item{statPredYv.max}{with the max algorithm, matrix indicating the statistical description of prediction categories for each observation on the validation datasets: number of predictions as an observation of the validation dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
\item{statPredYc.gravity}{with the gravity algorithm, matrix indicating the statistical description of prediction categories for each observation on the calibration datasets: number of predictions as an observation of the calibration dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
\item{statPredYv.gravity}{with the gravity algorithm, matrix indicating the statistical description of prediction categories for each observation on the validation datasets: number of predictions as an observation of the validation dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
\item{statPredYc.threshold}{with the threshold algorithm, matrix indicating the statistical description of prediction categories for each observation on the calibration datasets: number of predictions as an observation of the calibration dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
\item{statPredYv.threshold}{with the threshold algorithm, matrix indicating the statistical description of prediction categories for each observation on the validation datasets: number of predictions as an observation of the validation dataset, modal value, probability to be predicted with its standard deviation, 95\% confidence interval, quantiles 0.025 and 0.975, median value}
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)

Stone, M. (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the Royal Statistical Society B, 36(2), 111-147.
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\note{
at least 90 cross-validation repetitions may be recommended
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=plot_cvpred_mbplsda]{plot_cvpred_mbplsda}}
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
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
CVpred <- cvpred_mbplsda(modelembplsQ, nrepet = 30, threshold = 0.5, bloY = bloYobs, 
optdim = ncpopt, cpus = 1, algo = c("max"))

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
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
CVpred <- cvpred_mbplsda(modelembplsQ, nrepet = 90, threshold = 0.5, bloY = bloYobs, 
optdim = ncpopt, cpus = 1, algo = c("max"))
}

}

\keyword{multivariate}

