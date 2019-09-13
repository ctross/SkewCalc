#' Calculate the M index from the B index and relevant data.
#'
#' @param B The B index value estimated from a dataset.
#' @param R Total number of RS events in the sample from which B was estimated.
#' @param N Total number of individuals in the sample from which B was estimated.
#' @param t Assumption of the exposure time vector to be used in the correction term.
#' @param Samples Number of samples used to estimate the expected value of the correction term.
#' @return The M index.
#' @examples
#' set.seed(1)
#' B = 0.1 
#' R = 50
#' N = 10
#' M_index_from_B_index(B, R, N)

M_index_from_B_index = function(B, R, N, t=rep(1/N,N),Samples=1000) {
      E_Mraw = rep(NA,Samples)
    for(j in 1: Samples){
      t_hat = t/sum(t)
      E_Mraw[j] <- Mraw_index(rmultinom(1,R,t_hat),t)
      }
 M = Mraw_index_from_B_index(B,R,N) - mean(E_Mraw)
  return(M)
}

