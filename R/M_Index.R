#' M index of reproductive skew.
#'
#' @param r A vector of RS values.
#' @param t A vector of exposure times.
#' @param Samples Number of samples used to estimate the expected value of the correction term.
#' @return The M index.
#' @examples
#' set.seed(1) 
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index(RS, Age)

M_index = function(r, t, Samples=1000){
  if(min(t) <= 0){
  return(NA)
  }else{
    E_Mraw = rep(NA,Samples)
    for(j in 1: Samples){
      R = sum(r)
      t_hat = t/sum(t)
      E_Mraw[j] <- Mraw_index(rmultinom(1,R,t_hat),t)
      }
    M = Mraw_index(r,t) - mean(E_Mraw)
  return(M)
  }    
}
