
#' @title segfit -- LCE estimation of animal density from line transect survey data
#'
#' @description
#'  This function uses maximum likelihood estimation to estimate density (and some other related parameters) of animals from
#'  a mark-recapture line transect survey using two cameras, without any given recapture information. 
#'  The method is Latent Capture history Enumeration (LCE) and it is described in the paper 
#'  "A latent capture history model for digital aerial surveys" by D. L. Borchers, P. Nightingale, 
#'  B. C. Stevenson, and R. M. Fewster, to appear in the journal Biometrics.
#'  LCE uses a Markov model for animal availability (which in a marine survey can be interpreted 
#'  as a model of the diving cycle combined with a model of movements into or out of the detection strip).
#'  It uses a Gaussian movement model with a single speed parameter to model animal movement. 
#'  Since recaptures are unknown, LCE enumerates recapture scenarios. In each scenario where pairs of observations
#'  (that are sufficiently close) are 
#' 
#'  The package documentation has an example of using segfit.  
#'  
#' @param dat A list containing the following:
#' 
#' y1: Set of observations by the first observer, measured in distance units from the start of the transect.
#' 
#' y2: Set of observations by the second observer.
#' 
#' k: Lag (in seconds) between observer 1 and observer 2.
#' 
#' L: Length of the transect line, in the same distance units as y1 and y2. 
#' 
#' w: Half-width of the detection strip in distance units. 
#' 
#' b: Buffer half-width (distance from the centre line of the detection strip to the edge of the buffer), beyond which the method assumes that no animal can enter the detection strip between passage of the two observers. 
#' 
#' 
#' @param D.2D Starting point for the density estimate, expressed as a two-dimensional density. 
#' @param E1 Starting point for estimation of time in the surface or near-surface state during each dive cycle (in seconds).
#' @param Ec Dive cycle length in seconds.
#' @param sigmarate Starting point for estimation of speed of the animals (the parameter of the Gaussian movement model)
#' @param planespd Speed of the observers (in distance units per second)
#' @param p Vector of detection probabilities for the two observers (defaults to c(1,1))
#' @param sigma.mult  This parameter is no longer used (defaults to 5).
#' @param control.opt Passed through to the control.opt argument of function optim in package stats. 
#' @param method Optimisation method, defaults to "BFGS", passed to the optim function in package stats. 
#' @param estimate Vector of parameters to estimate, subset of c("D", "sigma", "E1", "mu_c") (D: density, sigma: animal speed, E1: time available for observation during dive cycle, mu_c: dive cycle length). Default value is c("D","sigma","E1"). 
#' @param set.parscale Set the parscale parameter of optim using the starting point values given (defaults to TRUE)
#' @param io Include in/out movement in the model, i.e. animals may become available or unavailable for detection between observers because of movement in/out of the detection strip (defaults to TRUE)
#' @param Dbound If estimating only D, Dbound may be used to provide upper and lower bounds for log(D). It is a list with elements $lower and $upper (default value NULL).
#' @param hessian If TRUE the Hessian matrix is calculated and returned. This can be used to estimate the asymptotic variance-covariance matrix of the parameters.
#' @param adj.mvt Adjust animal movement along the transect line to account for observer movement (e.g. if the animal moves in the opposite direction to the observers, then the elapsed time between the two observers passing the animal is less than k). If the flag is set to FALSE, the animal movement model assumes k seconds elapsed between the two observers passing the animal  (defaults to TRUE).
#' @param ft.normal If TRUE, uses normal to approximate Brownian hitting times, else uses exact expression for Brownian hitting times (defaults to FALSE).
#' @param cutstretch Factor to increase the maximum distance at which two observations may be considered to be observations of the same animal (defaults to 1).
#' @param krtest This parameter is no longer used. 


segfit=function(dat,D.2D,E1,Ec,sigmarate,planespd,p=c(1,1),sigma.mult=5,control.opt=NULL,
                method="BFGS",estimate=c("D","sigma","E1"),set.parscale=TRUE,
                io=TRUE,Dbound=NULL,hessian=TRUE,adj.mvt=TRUE,ft.normal=FALSE,cutstretch=1,krtest=FALSE) {
  
  require(car) # Needed for delta method se and cv of gamma
  
  k=dat$k
  s1 = dat$y1/planespd
  s2 = dat$y2/planespd
  tL = dat$L/planespd
  tw = dat$w/planespd
  tb = dat$b/planespd
  
  D.line.t=D.2D*2*dat$b*planespd # density in planespd along LINE units (1-dimensional)
  
  # Call C++ function to do the estimation
  
  fit=segfit.cpp(D=D.line.t,E1=E1,sigmarate=sigmarate,tw=tw,tb=tb,s1=s1,s2=s2,tL=tL,planespd=planespd,k=k,
                 p1=p[1],p2=p[2],Ec=Ec,hessian=hessian,control.opt=control.opt,method=method,estimate=estimate,
                 set.parscale-set.parscale,Dbound=Dbound,io=io,adj.mvt=adj.mvt,ft.normal=ft.normal,cutstretch=cutstretch)

  Dhat=fit$D/(2*dat$b*planespd)
  tau=fit$mu_c
  gamma=fit$E[1]/tau
  sigmarate=fit$sigmarate
  
  Dhat.se = gamma.se = sigmarate.se = NA
  Dhat.cv = gamma.cv = sigmarate.cv = NA
  Dhat.lcl = gamma.lcl = sigmarate.lcl = NA
  Dhat.ucl = gamma.ucl = sigmarate.ucl = NA
  vcv = NA
  if(hessian) {
    vcv=try(solve(fit$hessian),silent=TRUE)
    if(!inherits(vcv, "try-error")) {
      # Density
      Dhat.intest=logn.seci(log(fit$D),sqrt(vcv[1,1]))
      Dhat.se=Dhat.intest$se/(2*dat$b*planespd)
      Dhat.cv=Dhat.se/Dhat
      Dhat.lcl=Dhat.intest$lower/(2*dat$b*planespd)
      Dhat.ucl=Dhat.intest$upper/(2*dat$b*planespd)
      # gamma (need logit because it is logit(gamma), not gamma, that is actually estimated)
      gamma.intest=deltaMethod(c(x=logit(gamma)),"exp(x)/(1+exp(x))",vcov.=vcv[2,2])
      gamma.se=as.numeric(gamma.intest["SE"])
      gamma.cv=gamma.se/as.numeric(gamma.intest["Estimate"])
      gamma.lcl=as.numeric(gamma.intest[3])
      gamma.ucl=as.numeric(gamma.intest[4])
      # sigmarate
      sigmarate.intest=logn.seci(log(sigmarate),sqrt(vcv[3,3]))
      sigmarate.se=sigmarate.intest$se
      sigmarate.cv=sigmarate.se/sigmarate
      sigmarate.lcl=sigmarate.intest$lower
      sigmarate.ucl=sigmarate.intest$upper
    } 
  }
  
  D=c(est=Dhat,se=Dhat.se,cv=Dhat.cv,lcl=Dhat.lcl,ucl=Dhat.ucl)
  gamma=c(est=gamma,se=gamma.se,cv=gamma.cv,lcl=gamma.lcl,ucl=gamma.ucl)
  sigmarate=c(est=sigmarate,se=sigmarate.se,cv=sigmarate.cv,lcl=sigmarate.lcl,ucl=sigmarate.ucl)
  colnames(vcv) = row.names(vcv) = c("D","gamma","sigmarate")
  est = list(D=D, gamma=gamma, sigmarate=sigmarate, tau=tau,vcv=vcv)
  return(est)
}


segfit.cpp=function(D,E1,sigmarate,tw,tb,s1,s2,tL,planespd,k,p1,p2,Ec,hessian=FALSE,control.opt=NULL,method="BFGS",
                    estimate=c("D","E1","sigma"),set.parscale=TRUE,Dbound=NULL,io=FALSE,adj.mvt=FALSE,ft.normal=FALSE,
                    cutstretch=1)
  #-------------------------------------------------------------------------------
  # Just calls optim to maximise function rcpp_compute_likelihood() with respect 
  # to parameter D and g12 and/or sigma (suitably transformed to respect boundaries).
  #
  # This version segmentizes line for likelihood calculation.
  #
  # Inputs:
  # -------
# D         : density (Poisson process intensity parameter),
# E1        : Expected time available (in state 1) in seconds.
# sigmarate : std dev of normal governing animal movement in 1 sec. 
# dmax.t : maximum distance real duplicates could be apart (in plane seconds)
# s1   : locations (in plane seconds) of observer 1 detections
# s2   : locations (in plane seconds) of observer 2 detections
# tL    : total line length (in plane seconds)
# k    : lag (in seconds) between observer 1 and observer 2
# p1  : Bernoullis success probability for plane 1 observation process given availability,
# p2   : Bernoullis success probability for plane 1 observation process given availability,
# Ectrue : Expected value of dive cycle  (in seconds) = 1/g12+1/g21
# hessian: If TRUE, returns hessian
# control: optimization control (see help for optim())
# method: optimization method (see help for optim())
# estimate: character vector giving names of parameters to be estimated: 
#           subsets of c("D","E1","sigma","mu_c") are valid
# cutstretch: a factor by which to increas max dist that things could be considered pairs
# set.parscale: if TRUE sets parscale element of control list for optim
# Dbound: list with elements $lower and $upper, specifying bounds for log(D) if estimating only D
# io: if TRUE, accommodates in-out movement, else assumes no horizontal movement
# tw: strip half-width, in observer-seconds (only needed if io==TRUE)
#
# Outputs:
# -------
# List est with elements as per optim outputs, plus:
#   $D    : Estimated animal density
#   $G    : Estimated Markov transition matrix for availability process 
#           (g12 ie element (1,2) of this matrix.)
#   $E    : Expected dive time cycle length
#   $sigma: Estimated std dev of normal governing animal movement
#-------------------------------------------------------------------------------
#
#    NEW PARAMETERIZATION:
#    theta_1 = log of density
#    q_11 = Element 1,1 of the transition matrix.   q_11=-1/Ea  where Ea is the expected time available.
#    theta_3 = log of animal speed (sigma)
#    beta_2 = negative logit of gamma_21
{
  # set stuff up for returning:
  est=list(D=D,sigmarate=sigmarate,E=c(E1,Ec-E1),mu_c=Ec)
  # find segmentation cuts:
  dmax.t = tb-tw # max dist animal can move (in plane seconds)
  cuts=sort(unique(c(0,tL,segmentize(s1,s2,dmax.t*cutstretch))))
#  cuts=sort(unique(c(0,tL,segmentize(s1,s2,dmax.t))))
  
  if(length(estimate)==1 && estimate[1]!="D") stop("If only estimating one parameter, it must be D.")
  
  # build the theta vector. 
  theta=c()
  if(is.element("D",estimate)) {
    theta=append(theta, log(D))
  }
  if(is.element("E1",estimate)) {
#    theta=append(theta, log(E1))
    theta=append(theta, logit(E1/Ec))
  }
  if(is.element("sigma",estimate)) {
    theta=append(theta, log(sigmarate))
  }
  if(is.element("mu_c",estimate)) {
    theta=append(theta, log(Ec))
  }
  
  # set up default values 
  theta_1=log(D)
#  theta_2=log(E1)
  theta_2=logit(E1/Ec)
  theta_3=log(sigmarate)
  theta_4=log(Ec)
  
  if(set.parscale) control.opt$parscale=theta
  
  if(length(estimate)==1) {
    optout=optim(par=theta,fn=segnegllik.cpp.mix.io,s1=s1,s2=s2,dmax.t=dmax.t,p1=p1,p2=p2,k=k,planespd=planespd,theta_1=theta_1,
                 theta_2=theta_2,theta_3=theta_3,theta_4=theta_4,cuts=cuts,estimate=estimate,tw=tw,
                 control=control.opt,method="Brent",lower=log(D)-1,upper=log(D)+1,hessian=hessian,adj.mvt=adj.mvt,io=io,
                 ft.normal=ft.normal)
  } else {
    optout=optim(par=theta,fn=segnegllik.cpp.mix.io,s1=s1,s2=s2,dmax.t=dmax.t,p1=p1,p2=p2,k=k,planespd=planespd,theta_1=theta_1,
                 theta_2=theta_2,theta_3=theta_3,theta_4=theta_4,cuts=cuts,estimate=estimate,tw=tw,
                 control=control.opt,method=method,hessian=hessian,adj.mvt=adj.mvt,io=io,ft.normal=ft.normal)
  }
  
  # convert parameters to natural scale and save:
  
  # Get Ec first because need it to calculate E1 below
  if(is.element("mu_c",estimate)) {
    Ec=exp(optout$par[length(optout$par)]) # Ec is always last parameter
  } else {
    Ec = exp(theta_4)
  }
  
  if(is.element("D",estimate)) {
    if(io) {
      est$D=exp(optout$par[1])
    } else {
      est$D=exp(optout$par[1])
    }
    optout$par=optout$par[-1]
  }
  if(is.element("E1",estimate)) {
#    est$E[1]=exp(optout$par[1])
    est$E[1]=inv.logit(optout$par[1])*Ec
    optout$par=optout$par[-1]
  }
  if(is.element("sigma",estimate)) {
    est$sigmarate=exp(optout$par[1])
    optout$par=optout$par[-1]
  }
  if(is.element("mu_c",estimate)) {
    est$mu_c=exp(optout$par[1])
    optout$par=optout$par[-1]
  }
  
  est$E[2]=est$mu_c-est$E[1]
  
  # log-likelihood: 
  est$loglik=-optout$value
  est$AIC=-2*est$loglik+2*length(optout$par)
  est$hessian=NULL
  if(hessian) est$hessian=optout$hessian
  return(est)
}


segnegllik.cpp.mix.io=function(theta,s1,s2,dmax.t,p1,p2,k,planespd,adj.mvt=FALSE,
                               theta_1,theta_2,theta_3,theta_4,cuts,estimate,tw=0,io=TRUE,ft.normal=FALSE) {
  #-----------------------------------------------------------------------------------------
  # This is segnegllik.cpp.mix modified to deal with in-out movement
  # tw is strip half-width in observer-seconds
  # All of theta[1]=log(D), theta[2]=log(E1) and theta[3]=log(sigmarate) 
  # and theta[4]=log(Ec) are treated as parameters.
  # Unpack and compute likelihood given the parameters in 'estimate'
  #-----------------------------------------------------------------------------------------
  
  E1error=sigmaerror=FALSE # error trap for impossible/erroneous parameters
  
  theta_new=c()
  
  if(is.element("D",estimate)) {
    theta_new=append(theta_new, theta[1])
    theta=theta[-1]
  }
  else {
    theta_new=append(theta_new, theta_1)
  }
  
  if(is.element("E1",estimate)) {
    theta_new=append(theta_new, theta[1])
    theta=theta[-1]
  }
  else {
    theta_new=append(theta_new, theta_2)
  }
  
  if(is.element("sigma",estimate)) {
    theta_new=append(theta_new, theta[1])
    theta=theta[-1]
  }
  else {
    theta_new=append(theta_new, theta_3)
  }
  
  if(is.element("mu_c",estimate)) {
    theta_new=append(theta_new, theta[1])
    theta=theta[-1]
  }
  else {
    theta_new=append(theta_new, theta_4)
  }

  # unpack parameters
  ##    E1 = exp(theta_new[2]) # theta_new[2] is log of expected time available
  Ec = exp(theta_new[4]) # theta_new[4] is log of expected dive cycle length
  E1 = inv.logit(theta_new[2])*Ec # theta_new[2] is logit of propotion of time available
  sigmarate = exp(theta_new[3]) # theta_new[3] is log of sigmarate
  
  halfw.dist = tw*planespd # convert from observer seconds to km
  # Error trap: if sigma more than 100 times dmax.t*planespd, warn and return bad negative log likelihood:
  if(sigmarate*sqrt(k)>dmax.t*planespd) {
    warning("sigma > dmax.t*planespd: returning -Inf log-likelihood")
    sigmaerror=TRUE
  }
  # Error trap: if E1 > Ec:
  if(E1==Ec) {
    warning("E1 = Ec: returning -Inf log-likelihood")
    E1error=TRUE
  }
  # Error trap: if E1 -> 0:
  if(log(E1)<= -50) {
    warning("E1 < exp(-50): returning -Inf log-likelihood")
    E1error=TRUE
  }
  # Error trap: if Ec -> 0:
  if(log(Ec)<= -50) {
    warning("Ec < exp(-50): returning -Inf log-likelihood")
    E1error=TRUE
  }
  
  
  if(sigmaerror | E1error) {
    #loglik = Inf # penalise likelihood with implausibly large sigma
    #loglik = -1000000000    #  For the BFGS with bounding box method.
    loglik = -Inf
  } else {
    # get capture history probabilities:
    p=c(p1, p2)
    ps = p.t(E1,Ec,p,sigmarate,k,dmax.t,planespd,halfw.dist,adj.mvt,io,ft.normal=ft.normal)
    p10.k = ps$ch10
    p01.k = ps$ch01
    p11.k = ps$ch11
    
    nseg=length(cuts)-1
    
    sL_total=max(cuts)-min(cuts)
    
    loglik=0
    
    for(i in 1:nseg) {
      ss1=s1[cuts[i]<s1 & s1<=cuts[i+1]]
      ss2=s2[cuts[i]<s2 & s2<=cuts[i+1]]
      if(length(ss1)>0 | length(ss2)>0) {
        sL=cuts[i+1]-cuts[i]
        #      print(paste("In segment ",i," number of observations ",length(ss1),"  ",length(ss2),sep=""))
        loglik=loglik+rcpp_compute_likelihood(ss1, ss2, dmax.t, theta_new[1], theta_new[2], theta_new[3], theta_new[4], p1, p2, k, sL, planespd, p10.k, p01.k, p11.k)
      }
    }
  }
  #cat("Negative log likelihood:", -loglik, "\n")
  return(-loglik)
}


test.cpp.likelihood=function(s1, s2, dmax.time, l, l2, see.1, see.2, E1, E2, Ec, D, sigma, planespd, p1, p2, tL, dups, n1, n2, n11) {
  #  The C++ likelihood (given the true recaptures) and the known-recaptures likelihood should be the same
  #  UNLESS an animal was seen twice and moved more than dmax.time, in which case it is assumed by the C++ likelihood not to be a recapture.
  simt=l2-l  # Forward movement between planes. 
  
  varstate=c()
  for(i in 1:length(l)) {
    for(j in 1:length(l2)) {
      if(see.1[i]==1 && see.2[j]==1 && abs(l[i]-l2[j])<=dmax.time) {
        # l[i] and l2[j] are close enough to be potentially paired, and they have both been seen. 
        varstate=append(varstate, (i==j)*1)   # iff i==j they are the same animal.
      }
    }
  }
  
  
  # check varstate gives same number of (seen by 1), (seen by 2), (seen by both).
  
  ## check if one or other falls over when no recaptures. 
  
  q11=-1/(E1)    # Convert theta2 (log of Ea, expected time available) into q11
  q22=-1/(E2)    # Convert theta2 and theta4 (logs of Ea and expected dive cycle length) to q22
  p=c(p1, p2)
  Qmat=matrix(c(q11,-q22,-q11,q22),nrow=2)
  stdbn=c(Qmat[2,1],Qmat[1,2])/(Qmat[1,2]+Qmat[2,1])
  p01.k = p.omega.t(k,stdbn,p1,p2,Qmat,omega=01)
  p10.k = p.omega.t(k,stdbn,p1,p2,Qmat,omega=10)
  p11.k = p.omega.t(k,stdbn,p1,p2,Qmat,omega=11)
  
  ans1=calculate_Ljk(s1, s2, dmax.time, log(D), log(E1), log(sigma), log(Ec), p1, p2, k, tL, planespd, p10.k, p01.k, p11.k, varstate)
  
  n.omega=c(n1-n11, n2-n11, n11)
  ans2=-negllik(c(log(D), qlogis(E1/Ec), log(sigma/planespd)),k,n.omega,simt[dups],c(p1,p2),Ec,dmax.time,tL)
  
  if(sum(varstate)!=n11) print(paste("Warning: in test.cpp.likelihood sum of varstate is not same as m: ",sum(varstate),", ",n11,sep=""))
  # ans1 and ans2 should be the same apart from rounding error. 
  if(abs(ans1-ans2)>0.1 && sum(varstate)==n11) {
    print(paste("Warning: in test.cpp.likelihood C++ likelihood ",ans1," and known-recaptures likelihood ",ans2," do not match.",sep=""))
    #print(paste("sum(varstate):",sum(varstate)," varstate:",cat(varstate,sep=","),sep=""))
  }
}
