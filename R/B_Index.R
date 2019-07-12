#' B index a la Nonacs.
#'
#' @param r A vector of RS values.
#' @param t A vector of exposure times.
#' @return The B index.
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' B_index(RS, Age)

B_index = function(r, t){   
  if(min(t) <= 0){              
  return(NA)
  } 
  else{                      
  T = sum(t)
	Nbar = T / max(t)
	R = sum(r)
	if(R > 0){
	B = sum((r / R - t / T)^2) - (1 / R) * (1 - 1 / Nbar)
	} else {
	B = 0
	}
	return(B)
  }
}


