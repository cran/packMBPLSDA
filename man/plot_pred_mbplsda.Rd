\name{plot_pred_mbplsda}
\alias{plot_pred_mbplsda}

\title{
Plot the results of the fonction pred_mbplsda in a pdf file
}

\description{
Fonction to draw the results of the fonction pred_mbplsda (observed parameter values and predictions) in a pdf file
}

\usage{
plot_pred_mbplsda(obj, filename = "PlotPredMbplsda", propbestvar = 0.5)
}

\arguments{
  \item{obj}{
object type list containing the results of the fonction pred_mbplsda
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
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at 12emes Journees Scientifiques RFMF, Clermont-Ferrand, FRA(05-21-2019 - 05-23-2019).

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics, 15(10):134

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2020). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at Chimiometrie 2020, Liege, BEL(01-27-2020 - 01-29-2020).
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=mbplsda]{mbplsda}}
	\code{\link[=pred_mbplsda]{pred_mbplsda}}
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
bloYobs <- 2
ncpopt <- 1
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)
predictions <- pred_mbplsda(modelembplsQ, optdim = ncpopt, threshold = 0.5, 
bloY=bloYobs, algo = c("max", "gravity", "threshold"))
plot_pred_mbplsda(predictions,"plotPred_nf1", propbestvar=0.20)
}
}

\keyword{multivariate}

