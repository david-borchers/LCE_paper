\name{twoplane-package}
\alias{twoplane-package}
\alias{twoplane}
\docType{package}
\title{
LCE estimation for two-camera survey with unknown recaptures
}
\description{
Uses Latent Capture history Estimator (LCE) to estimate density from mark-recapture line transect survey using two cameras, without any recapture (duplicate) identification. Does this on the basis of the locations of detections by each observer, an animal movement model and a Markov model for the animal availability process.
}
\details{
\tabular{ll}{
Package: \tab twoplane\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-08-03\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
David Borchers and Peter Nightingale

Maintainer: David Borchers <dlb@st-andrews.ac.uk>
}
\references{
Borchers, D.L., Nightingale, P., Stevenson, B.C. and Fewster, R.M. (submitted) A latent capture history model for digital aerial surveys. 
}
\keyword{ package 
}
\seealso{
~~ Optional links to other man pages, e.g. \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(twoplane)
require(HiddenMarkov)
require(functional)
require(Rcpp)
require(boot)
require(expm)
require(Matrix)
require(mvtnorm)
require(car)
require(palm)


#  Get the porpoise data from package palm
data("porpoise")

# Set some constants and starting values
nm2km=1.852 # multiplier to convert nautical miles to kilometres
tau = 110 # dive cycle length
gamma = 86/110 # starting point for gamma estimate
kappa = gamma*tau # starting point for kappa estimate
D = 1.24 # starting point for density estimate
planeknots=100 # observer speed in knots
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
sigma=0.15*2/sqrt(248)   #  convert sigma from CCR to that for LCE
b = 2 # buffer half-width
w = 0.125 # striop half-width
sigma.mult = (b-w)/sigma # number of sigmas that b is greater than w


# Convert data to format required by segfit:
sdat = Palm2mleData(porpoise.data$points,porpoise.data$cameras,porpoise.data$d,porpoise.data$l,porpoise.data$w,b)

# Set optim() control parameters
control.opt=list(trace=0,maxit=1000)

# Specify which parameters to estimate:
# D = density,
# sigma = Brownian motion parameter
# E1 = expected time in state 1 (near-surface state)
estimate=c("D","sigma","E1") # parameters to estimate

# Get the MLE
mlefit<-segfit(sdat,D.2D=D,E1=kappa,Ec=tau,sigmarate=sigma,planespd=planespd,p=c(1,1),
               sigma.mult=sigma.mult,control.opt=control.opt,method="BFGS",estimate=estimate,
               set.parscale=TRUE,io=TRUE,Dbound=NULL,hessian=TRUE,adj.mvt=TRUE,ft.normal=FALSE)

# Look at the MLE
mlefit

# Look at the asymptotic estiamte of MLE parameter correlation
cov2cor(mlefit$vcv)

# Convert sigma to the sort of sigma that CCR method uses
sigma.palm = sigmarate2sigmapalm(mlefit$sigmarate["est"],248)

# Calculate mean velocity associated with sigma estimate
getspeed(mlefit$sigmarate["est"],248)*1000

}
