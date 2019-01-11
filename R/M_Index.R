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
  M = Mraw_index(r, t)*sqrt(mean(r))
  return(M)
  }
}
