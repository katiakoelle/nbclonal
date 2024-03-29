% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{Px}
\alias{Px}
\title{Probability of Extinction}
\usage{
Px(n, R0, mu, maxMuGen, maxFS)
}
\arguments{
\item{n}{Initial population size, a positive integer}

\item{R0}{Reproduction number, a positive number}

\item{mu}{Mutation rate, a positive number}

\item{maxMuGen}{Maximum number of mutant lineages being calculated}

\item{maxFS}{Maximum number of final sizes being calculated}
}
\value{
A number representing the probability
}
\description{
This function calculates the probability of the whole lineage goes extinction
through multiplying the probability that both the wild-type and the mutant
lineages goes extinction.
}
\examples{
Px(2, 1.2, 0.4, 60, 60)
Px(3, 3.5, 1.6, 120, 130)
}
