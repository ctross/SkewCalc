#' Calculate the M index from the B index and relevant data.
#'
#' @param B The B index value estimated from a dataset.
#' @param R Total number of RS events in the sample from which B was estimated.
#' @param N Total number of individuals in the sample from which B was estimated.
#' @param t Assumption of the exposure time vector used in the correction term.
#' @param Samples Samples used in estimating the correction term.
#' @return The M index.
#' @examples
#' B = 0.1 
#' R = 50
#' N = 10
#' M_index_from_B_index(B, R, N)

M_index_from_B_index = function(B, R, N, t=rep(1/N,N),Samples=1000) {
      E_Mraw_sq = rep(NA,Samples)
    for(j in 1: Samples){
      t_hat = t/sum(t)
      E_Mraw_sq[j] <- Mraw_index(rmultinom(1,R,t_hat),t)^2
      }
 M = Mraw_index_from_B_index(B,R,N)^2 - mean(E_Mraw_sq)
  return(M)
}

