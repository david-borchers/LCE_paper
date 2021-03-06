% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoplanesim.r
\name{dosim}
\alias{dosim}
\title{Simulates 2-plane data & fits LCE and CCR models to this}
\usage{
dosim(
  D.2D,
  L,
  w,
  b,
  sigmarate,
  k,
  planespd,
  kappa,
  tau,
  p = c(1, 1),
  movement = list(forward = TRUE, sideways = TRUE),
  fix.N = TRUE,
  En = NULL,
  Nsim = 100,
  writeout = TRUE,
  seed = 1,
  simethod = "MLE",
  control.opt = control.opt,
  hessian = TRUE,
  adj.mvt = FALSE,
  ft.normal = FALSE,
  sim.ft.normal = FALSE,
  progbar = TRUE,
  fn.append = NULL
)
}
\arguments{
\item{D.2D}{animal density (number per sq km)}

\item{L}{transect length (km)}

\item{w}{searched strip width (km)}

\item{b}{buffered strip width (km)}

\item{sigmarate}{the Brownian movement rate (sigma) parameter}

\item{k}{the lag between observers (seconds)}

\item{planespd}{Observer speed in km/sec}

\item{kappa}{mean time animals on surface (seconds)}

\item{tau}{mean dive cycle leength (seconds)}

\item{p}{vector of each observers' detection prob for available animals}

\item{movement}{list with components \code{forward} and \code{sideways}, which should be
TRUE or FALSE depending on what sort of movement you want to simulate with. 
Note: Estimation does not currently work when simulating with \code{forward=FALSE}.}

\item{fix.N}{Boolean variable, which if TRUE simulates with the same abundance, equal
to the expected, abundance, each simulation rep.}

\item{En}{Ignored if NULL, else the expected number of detections by each observer. 
Function sets transect length, \code{L}, to achieve this}

\item{Nsim}{Number of simulations to do}

\item{writeout}{Boolean, which if TRUE, causes results to be written to file and
the name of the file (which encodes the simulation parameters) to be returned
invisibly.}

\item{seed}{Integer seed which can be used for repoducible results}

\item{control.opt}{passed to \code{optim}}

\item{hessian}{Boolean controlling whether or not to return Hessian 
(and interval estimates for MLE method).}

\item{adj.mvt}{If TRUE, MLE method takes account of randomness in lag between
encounters; if FALSE, it sets these all to the lag, \code{k}.}

\item{ft.normal}{if TRUE, approximates Brownian motion hitting time with a 
normal distribution when estimation with MLE method, else uses exact expression
for Brownian hitting time}

\item{sim.ft.normal}{if TRUE, approximates Brownian motion hitting time with a 
normal distribution when simulating, else uses exact expression
for Brownian hitting time in simulation. 
NOTE: Exact simulator seems to have too low variance; with fast observers 
relative to animals, normal approx is very good (see function \code{dft} and
code in file \code{ft_plot.R}).}

\item{progbar}{If TRUE, puts a progress bar on screen so you can monitor simulation
progress.}

\item{fn.append}{A text string to append to the default output file name.}

\item{simthod}{Character; either `MLE` (for max likelihood simulator) or `Palm`
for CCR simulator from package \code{palm}.}
}
\description{
Simulates 2-plane data & fits LCE and CCR models to this
}
