#' Estimate the elasticity of RS on age.
#'
#' @param r A vector of RS values.
#' @param t A vector of ages at death or last census.
#' @param t A vector of ages at first census. Defaults to zero.
#' @return MLE of the elasticity of RS on Age, using Poisson regression.
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' elast(RS, Age)

elast = function(r,t,t0=0){
    f = function(x){
       d = exp(x[1])*(t^exp(x[2]) - t0^exp(x[2]))
       p = rep(NA,length(r))
       for(i in 1:length(r))
       p[i] = dpois(r[i], d[i], log = TRUE)
       return(-sum(p))
      }
   return(exp(optim( c(0.03,0.5), f )$par[2]))
}
