% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{pmf.P_clonal_round}
\alias{pmf.P_clonal_round}
\title{Round Probability Distribution of Getting Clonal Mutation}
\usage{
pmf.P_clonal_round(n, R0, mu, maxMuGen, maxFS, maxClonal)
}
\arguments{
\item{n}{Initial population size, a positive integer}

\item{R0}{Reproduction number, a positive number}

\item{mu}{Mutation rate, a positive number}

\item{maxMuGen}{Maximum number of mutant lineages being calculated}

\item{maxFS}{Maximum number of final sizes being calculated}

\item{maxClonal}{Maximum number of clonal mutations being calculated}
}
\value{
A data frame of clonal mutations and corresponding probabilities
}
\description{
This function calculates the round probability distribution of getting a
certain number of clonal mutations. Notice that this is not the exact
probability.
}
\examples{
pmf.P_clonal_round(2, 1.2, 0.4, 60, 60, 5)
pmf.P_clonal_round(3, 3.5, 1.6, 120, 130, 5)
}
