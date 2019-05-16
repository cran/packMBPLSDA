\name{testdim_mbplsda}
\alias{testdim_mbplsda}

\title{
Test of number of components by two-fold cross-validation for a multi-block partial least squares discriminant model
}

\description{
Function to perform a two-fold cross-validation in order to select the optimal number of dimensions of a multi-block partial least squares discriminant model, according to the classification error rate or to the area under ROC curve
}

\usage{
testdim_mbplsda(object, nrepet = 100, algo = c("max", "gravity", "threshold"),
threshold = 0.5, bloY, outputs = c("ER", "ConfMat", "AUC"), cpus = 1)
}

\arguments{
  \item{object}{
an object created by mbplsda_nfX
}
  \item{nrepet}{
integer indicating the number of repetitions
}
  \item{algo}{
character vector indicating the method(s) of prediction to use (see details)
}
  \item{threshold}{
numeric indicating the threshold, between 0 and 1, to consider the categories are predicted with the threshold prediction method.
}
  \item{bloY}{
integer vector indicating the number of categories per variable of the Y-block.
}
  \item{outputs}{
character vector indicating the wanted outputs (see details)
}
  \item{cpus}{
integer indicating the number of cpus to use when running the code in parallel
}

}

\details{
Three different algorithms are available to predict the categories of observations.
In the max, and respectively the threshold algorithms, numeric values are calculated from the matrix of explanatory variables and the regression coefficients. Then, the predicted categorie for each variable of the Y-block is the one which corresponds to the higher predicted value, respectively to the values higher than the indicated threshold.
In the gravity algorithm, predicted scores of the observations on the components are calculated. Then, each observation is assigned to the observed category of which it is closest to the barycentre in the component space.

Available outputs are Error Rates (ER), Confusion Matrix (ConfMat), Aera Under Curve (AUC).
}

\value{
	\item{TRUEnrepet}{number of repetitions}
	\item{TruePosC.max, .gravity, .threshold}{statistical description of percentages of true positive observations per category, evaluated on the calibration dataset, with the different algorithms (TPcM for "max", TPcG for "gravity", TPcT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{TruePosV.max, .gravity, .threshold}{statistical description of percentages of true positive observations per category, evaluated on the validation dataset, with the different algorithms (TPvM for "max", TPvG for "gravity", TPvT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{TrueNegC.max, .gravity, .threshold}{statistical description of percentages of true negative observations per category, evaluated on the calibration dataset, with the different algorithms (TNcM for "max", TNcG for "gravity", TNcT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{TrueNegV.max, .gravity, .threshold}{statistical description of percentages of true negative observations per category, evaluated on the validation dataset, with the different algorithms (TNvM for "max", TNvG for "gravity", TNvT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{FalsePosC.max, .gravity, .threshold}{statistical description of percentages of false positive observations per category, evaluated on the calibration dataset, with the different algorithms (FPcM for "max", FPcG for "gravity", FPcT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{FalsePosV.max, .gravity, .threshold}{statistical description of percentages of false positive observations per category, evaluated on the validation dataset, with the different algorithms (FPvM for "max", FPvG for "gravity", FPvT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{FalseNegC.max, .gravity, .threshold}{statistical description of percentages of false negative observations per category, evaluated on the calibration dataset, with the different algorithms (FNcM for "max", FNcG for "gravity", FNcT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{FalseNegV.max, .gravity, .threshold}{statistical description of percentages of false negative observations per category, evaluated on the validation dataset, with the different algorithms (FNvM for "max", FNvG for "gravity", FNvT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{ErrorRateC.max, .gravity, .threshold}{statistical description of prediction error rates per category, evaluated on the calibration dataset, with the different algorithms (ERcM for "max", ERcG for "gravity", ERcT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{ErrorRateV.max, .gravity, .threshold}{statistical description of prediction error rates per category, evaluated on the validation dataset, with the different algorithms (ERvM for "max", ERvG for "gravity", ERvT for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{ErrorRateCglobal.max, .gravity, .threshold}{statistical description of global prediction error rates, evaluated on the calibration dataset, with the different algorithms (ERcM.global for "max", ERcG.global for "gravity", ERcT.global for "threshold"), for a number of components ranging from 1 to its maximum value}
	\item{ErrorRateVglobal.max, .gravity, .threshold}{statistical description of global prediction error rates, evaluated on the validation dataset, with the different algorithms (ERvM.global for "max", ERvG.global for "gravity", ERvT.global for "threshold"), for a number of components ranging from 1 to its maximum value}
	
	\item{AUCc}{statistical description of aera under ROC curve values per category, evaluated on the calibration dataset, if all Y-block variables are binary, for a number of components ranging from 1 to its maximum value}
	\item{AUCv}{statistical description of aera under ROC curve values per category, evaluated on the validation dataset, if all Y-block variables are binary, for a number of components ranging from 1 to its maximum value}
	\item{AUCc.global}{statistical description of global aera under ROC curve values, evaluated on the calibration dataset, if all Y-block variables are binary, for a number of components ranging from 1 to its maximum value}
	\item{AUCv.global}{statistical description of global aera under ROC curve values, evaluated on the validation dataset, if all Y-block variables are binary, for a number of components ranging from 1 to its maximum value}
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
at least 30 cross-validation repetitions may be recommended
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=plot_testdim_mbplsda]{plot_testdim_mbplsda}}
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
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 3)
resdim <- testdim_mbplsda(object = modelembplsQ, nrepet = 30, threshold = 0.5, 
bloY = bloYobs, cpus = 1, algo = c("max"), outputs = c("ER"))
}

\keyword{multivariate}

