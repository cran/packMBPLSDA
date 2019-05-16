\name{disjunctive}
\alias{disjunctive}

\title{
Disjunctive table
}

\description{
Function to transform a boolean matrix in a disjunctive table
}

\usage{
disjunctive(y)
}

\arguments{
  \item{y}{boolean matrix indicating observations categories}
}
\details{
no details are needed
}

\value{
\item{ydisj}{disjunctive table}
}

\references{
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. 12emes Journees Scientifiques RFMF (05-21-2019 - 05-23-2019) Clermont-Ferrand(FRA). In: Livre des resumes RFMF 2019 (in press)

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics (submitted)
}

\author{
Marion Brandolini-Bunlon (<marion.brandolini-bunlon@inra.fr>) and Stephanie Bougeard (<stephanie.bougeard@anses.fr>)
}

\seealso{
	\code{\link[=packMBPLSDA-package]{packMBPLSDA-package}}
}

\examples{
data(status)
disjonctif <- (disjunctive(status))
}

\keyword{datagen}