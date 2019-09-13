#' M index of reproductive skew using Stan to conduct Bayesian bootstrapping
#'
#' @param r A vector of RS values over time period of observation (i.e., between times t0 and t).
#' @param t A vector of ages at death or last census.
#' @param t0 A vector of ages at first census. Defaults to zero.
#' @param samples Number of MCMC samples per chain.
#' @param warmup Number of warmup iterations per chain.
#' @param chains Number of chains.
#' @param adapt_delta A tuning parameter in Stan. See "rstan" package for more details.
#' @param max_treedepth A tuning parameter in Stan. See "rstan" package for more details.
#' @param refresh How often to print updates on MCMC progress.
#' @return A Stan object "StanResults" and a data object "model_dat" are saved to the workspace. These can be used for plots and calculations. Results are also printed: 
#' @return 1) M. The posterior distribution of M, not accounting for diminishing RS returns to age. 
#' @return 2) M_age. The posterior distribution of M, accounting for diminishing RS returns to age. 
#' @return 3) M_raw. The posterior distribution of M_raw, not accounting for diminishing RS returns to age. 
#' @return 4) M_raw_age. The posterior distribution of M_raw, accounting for diminishing RS returns to age. 
#' @return 5) gamma The estimated elasticity of RS on exposure time. 
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index_stan(RS, Age)


M_index_stan = function(r, t, t0=FALSE, samples=2000, warmup=1000, chains=2, adapt_delta=0.9, max_treedepth=12, refresh=1) {
  if(min(t)<=0){
   return(NA)
   }else{
   if(length(t0)==1)
   t0 = rep(0,length(r))

   model_dat<-list(
    N=length(r),
    r=r,
    t=t,
    t0=t0
    )
    
    model_dat<<-model_dat
    
StanResults <<- stan(model_code=model_code, data=model_dat, thin=1, iter=samples, 
                    warmup=warmup, chains=chains, refresh=refresh,
                    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth))

print(StanResults,pars=c("M", "M_age","M_raw", "M_raw_age", "gamma"))
}
}
