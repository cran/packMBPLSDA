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
Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at 12emes Journees Scientifiques RFMF, Clermont-Ferrand, FRA(05-21-2019 - 05-23-2019).

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2019). Multi-block PLS discriminant analysis for the joint analysis of metabolomic and epidemiological data. Metabolomics, 15(10):134

Brandolini-Bunlon, M., Petera, M., Gaudreau, P., Comte, B., Bougeard, S., Pujos-Guillot, E.(2020). A new tool for multi-block PLS discriminant analysis of metabolomic data: application to systems epidemiology. Presented at Chimiometrie 2020, Liege, BEL(01-27-2020 - 01-29-2020).
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