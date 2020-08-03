#' @title 95\% Std deviation and confidence interval for lognormal
#'
#' @description
#'  Calculates the standard deviation and 95% confidence interval for a lognormal
#'  random variable, given the normal mean and standard deviation.
#'  
#' @param mu Normal mean
#' @param sd Normal std dev.
#' @examples 
#' sigma.beta=0.24
#' mu.beta= -1.6
#' logn.seci(mu.beta,sigma.beta)
#' 
logn.seci=function(mu,sd) {
  se = sqrt((exp(sd^2)-1)*exp(2*mu + sd^2))
  bnd = exp(mu+c(-1,1)*1.96*sd)
  return(list(se=se,lower=bnd[1],upper=bnd[2]))
}