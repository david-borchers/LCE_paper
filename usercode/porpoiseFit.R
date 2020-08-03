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

# Look at numbers within segments:
# segmentplot(sdat,planespd)

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

