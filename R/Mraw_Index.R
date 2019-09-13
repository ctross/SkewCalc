#' Mraw index of reproductive skew.
#'
#' @param r A vector of RS values.
#' @param t A vector of exposure times.
#' @return The Mraw index.
#' @examples
#' set.seed(1) 
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' Mraw_index(RS, Age)

Mraw_index = function(r, t){
  if(min(t) <= 0){
  return(NA)
  }
  else{
	N = length(r)
	R = sum(r)
	T = sum(t)
	si = ((r / R) - (t / T))^2
	S = sum(si)
	Mraw = N * S
  return(Mraw)
  }
}
