\name{permut_mbplsda}
\alias{permut_mbplsda}

\title{
Permutation testing of a multi-block partial least squares discriminant model
}

\description{
Function to perform permutation testing with 2-fold cross-validation for multi-block partial least squares discriminant analysis, in order to evaluate model validity and predictivity
}

\usage{
permut_mbplsda(object, optdim, bloY, algo = c("max", "gravity", "threshold"), 
threshold = 0.5, nrepet = 100, npermut = 100, nbObsPermut = NULL, 
outputs = c("ER", "ConfMat", "AUC"), cpus = 1)
}

\arguments{
  \item{object}{
an object created by mbplsda_nfX
}
  \item{optdim}{
integer indicating the (optimal) number of components of the multi-block partial least squares discriminant model
}
  \item{bloY}{
integer vector indicating the number of categories per variable of the Y-block.
}
  \item{algo}{
character vector indicating the method(s) of prediction to use (see details)
}
  \item{threshold}{
numeric indicating the threshold, between 0 and 1, to consider the categories are predicted with the threshold prediction method.
}
  \item{nrepet}{
integer indicating the number of repetitions
}
  \item{npermut}{
integer indicating the number of Y-block with switching observations
}
  \item{nbObsPermut}{
integer indicating the number of switching observations in all the modified Y-blocks
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

If nbObsPermut is not NULL, t-test are performed to compare mean cross-validated overall prediction error rates (or aera under ROC curve) evaluated on permuted Y-blocks, with the cross-validated overall prediction error rate (or aera under ROC curve) evaluated on the original Y-block.

Available outputs are Error Rates (ER), Confusion Matrix (ConfMat), Aera Under Curve (AUC).
}

\value{
	\item{RV.YYpermut.values}{RV coefficient between Y-block and each Y-block with permuted values}
	\item{cor.YYpermut.values}{correlation coefficient between categories in the Y-block and each Y-block with permuted values}
	\item{prctGlob.Ychange.values}{overall percentage of modified values in each Y-block with permuted values}
	\item{prct.Ychange.values}{percentage per category of modified values in each Y-block with permuted values}
	\item{descrYperm}{statistical description of RV.YYpermut, cor.YYpermut, prctGlob.Ychange, prct.Ychange}
	
	\item{TruePosC.max, TruePosC.gravity, TruePosC.threshold}{statistical description of cross-validated percentages of true positive observations per category, evaluated on calibration datasets, with the different algorithms (TruePosC.max for "max", TruePosC.gravity for "gravity", TruePosC.threshold for "threshold"), for each Y-block with permuted values}
	\item{TruePosV.max, TruePosV.gravity, TruePosV.threshold}{statistical description of cross-validated percentages of true positive observations per category, evaluated on validation datasets, with the different algorithms (TruePosV.max for "max", TruePosV.gravity for "gravity", TruePosV.threshold for "threshold"), for each Y-block with permuted values}
	\item{TrueNegC.max, TrueNegC.gravity, TrueNegC.threshold}{statistical description of cross-validated percentages of true negative observations per category, evaluated on calibration datasets, with the different algorithms (TrueNegC.max for "max", TrueNegC.gravity for "gravity", TrueNegC.threshold for "threshold"), for each Y-block with permuted values}
	\item{TrueNegV.max, TrueNegV.gravity, TrueNegV.threshold}{statistical description of cross-validated percentages of true negative observations per category, evaluated on validation datasets, with the different algorithms (TrueNegV.max for "max", TrueNegV.gravity for "gravity", TrueNegV.threshold for "threshold"), for each Y-block with permuted values}
	\item{FalsePosC.max, FalsePosC.gravity, FalsePosC.threshold}{statistical description of cross-validated percentages of false positive observations per category, evaluated on calibration datasets, with the different algorithms (FalsePosC.max for "max", FalsePosC.gravity for "gravity", FalsePosC.threshold for "threshold"), for each Y-block with permuted values}
	\item{FalsePosV.max, FalsePosV.gravity, FalsePosV.threshold}{statistical description of cross-validated percentages of false positive observations per category, evaluated on validation datasets, with the different algorithms (FalsePosV.max for "max", FalsePosV.gravity for "gravity", FalsePosV.threshold for "threshold"), for each Y-block with permuted values}
	\item{FalseNegC.max, FalseNegC.gravity, FalseNegC.threshold}{statistical description of cross-validated percentages of false negative observations per category, evaluated on calibration datasets, with the different algorithms (FalseNegC.max for "max", FalseNegC.gravity for "gravity", FalseNegC.threshold for "threshold"), for each Y-block with permuted values}
	\item{FalseNegV.max, FalseNegV.gravity, FalseNegV.threshold}{statistical description of cross-validated percentages of false negative observations per category, evaluated on validation datasets, with the different algorithms (FalseNegV.max for "max", FalseNegV.gravity for "gravity", FalseNegV.threshold for "threshold"), for each Y-block with permuted values}
	\item{ErrorRateC.max, ErrorRateC.gravity, ErrorRateC.threshold}{statistical description of cross-validated prediction error rates per category, evaluated on calibration datasets, with the different algorithms (ErrorRateC.max for "max", ErrorRateC.gravity for "gravity", ErrorRateC.threshold for "threshold"), for each Y-block with permuted values}
	\item{ErrorRateV.max, ErrorRateV.gravity, ErrorRateV.threshold}{statistical description of cross-validated prediction error rates per category, evaluated on validation datasets, with the different algorithms (ErrorRateV.max for "max", ErrorRateV.gravity for "gravity", ErrorRateV.threshold for "threshold"), for each Y-block with permuted values}
	\item{ErrorRateCglobal.max, ErrorRateCglobal.gravity, ErrorRateCglobal.threshold}{statistical description of cross-validated overall prediction error rates, evaluated on calibration datasets, with the different algorithms (ErrorRateCglobal.max for "max", ErrorRateCglobal.gravity for "gravity", ErrorRateCglobal.threshold for "threshold"), for each Y-block with permuted values}
	\item{ErrorRateVglobal.max, ErrorRateVglobal.gravity, ErrorRateVglobal.threshold}{statistical description of cross-validated overall prediction error rates, evaluated on validation datasets, with the different algorithms (ErrorRateVglobal.max for "max", ErrorRateVglobal.gravity for "gravity", ErrorRateVglobal.threshold for "threshold"), for each Y-block with permuted values}
	
	\item{AUCc}{if all Y-block variables are binary, statistical description of cross-validated aera under ROC curve values per category, evaluated on the validation datasets, for each Y-block with permuted values}
	\item{AUCv}{if all Y-block variables are binary, statistical description of cross-validated aera under ROC curve values per category, evaluated on the validation datasets, for each Y-block with permuted values}
	\item{AUCc.global}{if all Y-block variables are binary, statistical description of cross-validated overall aera under ROC curve values, evaluated on the validation datasets, for each Y-block with permuted values}
	\item{AUCv.global}{if all Y-block variables are binary, statistical description of cross-validated overall aera under ROC curve values, evaluated on the validation datasets, for each Y-block with permuted values}

	\item{reg.GlobalRes_prctYchange}{results of linear regression of overall prediction error rates, and overall aera under ROC curve, onto percentages of modified values in Y-block}
	\item{ttestMeanERv}{if nbObsPermut is not NULL, results of the t-test comparing mean cross-validated overall prediction error rates (and eventually aera under ROC curve) evaluated on permuted Y-blocks, with the cross-validated overall prediction error rate (and eventually aera under ROC curve) evaluated on the original Y-block}
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

\note{
at least 30 cross-validation repetitions and 100 Y-block with switching observations may be recommended
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=plot_permut_mbplsda]{plot_permut_mbplsda}}
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
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
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", 
scannf = FALSE, nf = 1)
rtsPermut <- permut_mbplsda(modelembplsQ, nrepet = 30, npermut = 100, optdim = ncpopt, 
outputs = c("ER"), bloY = bloYobs, nbObsPermut = 10, cpus=1, algo = c("max"))
}
}

\keyword{multivariate}

