% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalMuFunctions.R
\name{P_wt_est}
\alias{P_wt_est}
\title{Wild-type Establishment Probability}
\usage{
P_wt_est(n, R0, mu)
}
\arguments{
\item{n}{Initial population size, a positive integer}

\item{R0}{Reproduction number, a positive number}

\item{mu}{Mutation rate, a positive number}
}
\value{
A number representing the probability
}
\description{
This function calculates the probability of establishment of wild-type
offspring of a lineage starting with certain parameters.
}
\examples{
P_wt_est(2, 1.2, 0.4)
P_wt_est(1, 3.5, 1.6)
}
