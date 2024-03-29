% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{pmf.P_clonal}
\alias{pmf.P_clonal}
\title{Probability Distribution of Getting Clonal Mutation Normalized}
\usage{
pmf.P_clonal(n, R0, mu, maxMuGen, maxFS, maxClonal)
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
This function calculates the probability distribution of getting a certain
number of clonal mutations. The resulting probabilities are normalized by
establishment, 1.e. 1-Px, and Px will not be contained in the final data
frame.
}
\examples{
pmf.P_clonal(2, 1.2, 0.4, 60, 60, 5)
pmf.P_clonal(3, 3.5, 1.6, 120, 130, 8)
}
