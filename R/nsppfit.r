twoplane.fit=function(sdat,tau,R,all=TRUE){
  pts=matrix(c(sdat$y1,sdat$y2),ncol=1)
  cameras=c(rep(1,length(sdat$y1)),rep(2,length(sdat$y2)))
  d=sdat$L
  w=sdat$w
  b=sdat$b
  l=sdat$k

  est = fit.twocamera(points=pts, cameras=cameras, d=d, w=w, b=b, l=l, tau=tau, R=R)
  
  params=coef(est)
  
  if(all) return(est) else return(params)
}

#' @title Palm2mleData -- Convert data set from Palm format to the format used by the LCE method. 
#' @param points Set of observations by both observers, measured in distance units from the start of the transect line. 
#' @param cameras Vector of integers 1 or 2. For each observation points[i], cameras[i] indicates which camera made the observation. 
#' @param d Length of the transect line, in the same distance units as points. 
#' @param l Lag (in seconds) between observer 1 and observer 2.
#' @param w Half-width of the detection strip in distance units. 
#' @param b Buffer half-width (distance from the centre line of the detection strip to the edge of the buffer), beyond which the method assumes that no animal can enter the detection strip between passage of the two observers. 
#' @return A list in the format required for the first argument of the segfit function.

Palm2mleData = function(points,cameras,d,l,w,b) {
  y1 = points[cameras==1]
  y2 = points[cameras==2]
  sdat = list(y1=y1,y2=y2,k=l,L=d,w=w,b=b)
  return(sdat)
}
