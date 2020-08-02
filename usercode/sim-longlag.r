
## Set simulation set number
i=1
## Specify number of simulation
Nsim=1
## Set seed for reprodicibility.
set.seed(i)

## Get the LCE and CCR estimation packages
library(twoplane) # LCE
library(palm) # CCR


nm2km=1.852 # multiplier to convert nautical miles to kilometres
planeknots=100 # observer speed in knots   
planespd=planeknots*nm2km/(60^2) # observer speed in km/sec

D = D.2D <- 1.05
## Time between cameras (seconds).
k = l <- 248
## Mean dive-cycle duration (seconds).
tau <- 110
## Mean duration of surface phase (seconds).
kappa <- 0.86*tau
sigma_palm = 0.15 # estimated Palm-type sigma (in km) from porpoise data, with lag 248 seconds
sigmarate = sigmapalm2sigmarate(sigma_palm,lag=248) # Brownian sigmarate consistent with sigma_palm
sigma = sigmarate*sqrt(k) # Standard error of Brownian movement after k seconds
animalspd = getspeed(sigmarate,248)*1000 # average speed of animals in m/s, after 248 seconds
planespd/(animalspd/1000) # How much faster plane is going than average animal speed
#speed2sigmarate(.95/1000,248) # Bownian motion parameter that gives the observed speed over 248 seconds of .95 m/s
sigma.mult=5 # multiplier of sigma, to set bound for maximum distance moved in k seconds

p=c(1,1) # Probability see, given available in searched strips

## Transect half-width.
halfw.dist = w <- 0.125
## Buffer distance
b <- w + sigma.mult*sigma

control.opt=list(trace=0,maxit=1000)
estimate=c("D","sigma","E1") # parameters to estimate
movement = list(forward=TRUE,sideways=TRUE)

# Here to simulate wiht Palm simulator:
# adj.mvt=TRUE makes LCE estimator allow time between encouners to depend on animal movement (not be fixed at k)
# ft.normal=FALSE makes LCE estimator use exact expression for Brownian hitting time when adj.mvt==TRUE
# palmvt$file contains names of files that contain simulation results after dosim has finished
progbar = FALSE # set to FALSE if you don't want to generate a bar to track progress of simulations
fn.append = paste("set", i, sep = "")

# Here to simulate wiht LCE simulator:
# sim.ft.normal=TRUE makes LCE simulator use normal approxx for Brownian hitting time 
# mlemvt$file contains names of files that contain simulation results after dosim has finished
mletmvt = system.time(mlemvt <- dosim(D.2D,L,w,b,sigmarate,k,planespd,kappa,tau,p=p,movement=movement,
                                      fix.N=FALSE,En=NULL,Nsim=Nsim,writeout=TRUE,seed=seed,simethod="MLE",
                                      control.opt=control.opt,adj.mvt=TRUE,ft.normal=FALSE,sim.ft.normal=TRUE,
                                      progbar=progbar,fn.append=fn.append))

mletmvt
