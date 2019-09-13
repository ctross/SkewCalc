#' Mraw index accounting for diminishing returns to age.
#'
#' @param r A vector of RS values.
#' @param t A vector of ages at death or last census.
#' @param t A vector of ages at first census. Defaults to zero.
#' @return The Mraw index accounting for diminishing returns to age.
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' Mraw_index_age(RS, Age)

Mraw_index_age = function(r, t, t0=0) {
  if(min(t-t0) <= 0){
   return(NA)
   }else{
   N = length(r)
   R = sum(r)
   beta = elast(r,t,t0) 
   tt = t^beta - t0^beta
   T = sum(tt)
   si = ((r/R)-(tt/T))^2
   S = sum(si)
   Mraw = N * S
   return(Mraw)
  }
}
