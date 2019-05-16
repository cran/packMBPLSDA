\name{ginv}
\alias{ginv}

\title{
generalized inverse of a matrix X
}

\description{
function to calculate the generalized inverse of a matrix X
}

\usage{
ginv(X, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{X}{
Matrix for which the generalized inverse is required
}
  \item{tol}{
A relative tolerance to detect zero singular values
}
}
