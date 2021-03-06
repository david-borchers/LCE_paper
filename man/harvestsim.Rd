% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoplanesim.r
\name{harvestsim}
\alias{harvestsim}
\title{Gets, and summarises simulation results}
\usage{
harvestsim(fn, badcut = 100)
}
\arguments{
\item{fn}{Name of the file with the output from \code{dosim}.
You should use the filename returned by \code{dosim};  see \code{dosim}
code for details.}
}
\description{
Gets, and summarises simulation results
}
\details{
Reads simulation results from \code{dosim}, using appropriate file name 
 (see details in code) and returns summary of bias, cv, coverage etc.
}
