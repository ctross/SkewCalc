#' M index accounting for diminishing returns to age.
#'
#' @param r A vector of RS values.
#' @param t A vector of ages at death or last census.
#' @param t A vector of ages at first census. Defaults to zero.
#' @return The M index accounting for diminishing returns to age.
#' @examples
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index_age(RS, Age)

M_index_age = function(r,t,t0=0) {
  if(min(t-t0) <= 0){
  return(NA)
   } else{
     E_Mraw_sq = rep(NA,Samples)
    for(j in 1: Samples){
      R = sum(r)
      t_hat = t/sum(t)
      E_Mraw_sq[j] <- Mraw_index_age(rmultinom(1,R,t_hat),t,t0)^2
      }
    M = Mraw_index_age(r,t,t0)^2 - mean(E_Mraw_sq)
  return(M)
  }
}
