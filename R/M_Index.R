#' M index of reproductive skew.
#'
#' @param r A vector of RS values.
#' @param t A vector of exposure times.
#' @return The M index.
#' @examples
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index(RS, Age)

M_index = function(r, t){
  if(min(t) <= 0){
  return(NA)
  }else{
    X = rmultinom(1,sum(r),t/sum(t))
    M = Mraw_index(r,t)^2 - Mraw_index(X,t)^2
  return(M)
  }    
}
