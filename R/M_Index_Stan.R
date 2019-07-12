#' M index of reproductive skew using Stan to simulate the population distribution of RS and age.
#'
#' @param r A vector of RS values.
#' @param t A vector of ages at death or last census.
#' @param t0 A vector of ages at first census. Defaults to zero.
#' @param Fast If TRUE then the model uses a faster model specification. If FALSE, then the model uses truncated distributions to more accuratly represent the data distribution.
#' @return A Stan object is saved to the workspace. This will be used for plots and calculations. The M and Mraw results have 4 columns each: 
#' @return 1) The point estimate on the sample data. 
#' @return 2) The point estimate on the sample data, but adjusting for diminising fitness returns to age.
#' @return 3) The posterior estimate using data simulated from the generative model. 
#' @return 4) The posterior estimate using data simulated from the generative model, but adjusting for diminising fitness returns to age.
#' @examples
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index_stan(RS, Age)

M_index_stan = function(r, t, t0=FALSE, Samples=2000, Warmup=1000, Chains=2, adapt_delta=0.9, max_treedepth=12, refresh=1) {
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
    
StanResults <<- stan(model_code=model_code, data=model_dat, thin=1, iter=Samples, 
                    warmup=Warmup, chains=Chains, refresh=refresh,
                    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth))

print(StanResults,pars=c("M", "M_age","M_raw", "M_raw_age", "gamma"))
}
}
