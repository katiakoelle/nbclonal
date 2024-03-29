% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{P_clonal}
\alias{P_clonal}
\title{Probability of Getting Clonal Mutation}
\usage{
P_clonal(n, m, R0, mu, maxMuGen, maxFS)
}
\arguments{
\item{n}{Initial population size, a positive integer}

\item{m}{Number of clonal mutation, a non-negative integer}

\item{R0}{Reproduction number, a positive number}

\item{mu}{Mutation rate, a positive number}

\item{maxMuGen}{Maximum number of mutant lineages being calculated}

\item{maxFS}{Maximum number of final sizes being calculated}
}
\value{
A number representing the probability
}
\description{
This function calculates the round probability of getting a certain number
of clonal mutations.
}
\examples{
P_clonal(2, 4, 1.2, 0.4, 60, 60)
P_clonal(3, 0, 3.5, 1.6, 120, 130)
}
