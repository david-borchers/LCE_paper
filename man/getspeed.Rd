% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma_speed.r
\name{getspeed}
\alias{getspeed}
\title{Average speed of Brownian motion particle}
\usage{
getspeed(sigmarate, lag = 1)
}
\arguments{
\item{sigmarate}{The Brownian motion with movement rate paramter}

\item{lag}{The time over which the average speed is required.}
}
\description{
Calculates the average speed E(displacement/time) of a particle doing 
 Brownian motion with movement rate paramter \code{sigmarate}, over a time
 period \code{lag}. 
 Uses the fact that the expected value of a Chi random variable is sqrt(2)/gamma(0.5).
}
