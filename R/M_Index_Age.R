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
    X = rmultinom(1,sum(r),t/sum(t))
    M = Mraw_index_age(r,t,t0)^2 - Mraw_index_age(X,t,t0)^2
  return(M)
  }
}
