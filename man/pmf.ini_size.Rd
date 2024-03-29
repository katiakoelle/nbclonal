% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{pmf.ini_size}
\alias{pmf.ini_size}
\title{Initial Population Size Distribution}
\usage{
pmf.ini_size(lambda, R0, maxIni)
}
\arguments{
\item{lambda}{Mean of Poisson distribution, a positive number}

\item{R0}{Reproduction number, a positive number}

\item{maxIni}{Maximum initial population size being calculated}
}
\value{
A data frame of N, distribution type and corresponding probability of getting that number of initial population size
}
\description{
This function calculates the poisson Distribution of initial population size
N given a mean of lambdaand also the adjusted distribution normalized by at
least one of the initial viral particles established.
}
\examples{
pmf.ini_size(2.1, 1.6, 10)
pmf.ini_size(0.8, 2.4, 6)
}
