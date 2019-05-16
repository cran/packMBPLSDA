\name{packMBPLSDA-package}
\alias{packMBPLSDA-package}
\alias{packMBPLSDA}
\docType{package}

\title{
\packageTitle{packMBPLSDA}
}

\description{
\packageDescription{packMBPLSDA}
}

\details{
\packageIndices{packMBPLSDA}
}

\author{

\packageAuthor{packMBPLSDA}

Maintainer: \packageMaintainer{packMBPLSDA}
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)
}

\keyword{package}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=testdim_mbplsda]{testdim_mbplsda}}
	\code{\link[=plot_testdim_mbplsda]{plot_testdim_mbplsda}}
	\code{\link[=permut_mbplsda]{permut_mbplsda}}
	\code{\link[=plot_permut_mbplsda]{plot_permut_mbplsda}}
	\code{\link[=pred_mbplsda]{pred_mbplsda}}
	\code{\link[=plot_pred_mbplsda]{plot_pred_mbplsda}}
	\code{\link[=cvpred_mbplsda]{cvpred_mbplsda}}
	\code{\link[=plot_cvpred_mbplsda]{plot_cvpred_mbplsda}}
	\code{\link[=boot_mbplsda]{boot_mbplsda}}
	\code{\link[=plot_boot_mbplsda]{plot_boot_mbplsda}}
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
