#' Calculate the Mraw index from the B index and relevant data.
#'
#' @param B The B index value estimated from a dataset.
#' @param R Total number of RS events in the sample from which B was estimated.
#' @param N Total number of individuals in the sample from which B was estimated.
#' @return The Mraw index.
#' @examples
#' B = 0.1 
#' R = 50
#' N = 10
#' Mraw_index_from_B_index(B, R, N)

Mraw_index_from_B_index = function(B,R,N) {
  M = B*N + (N-1)/R
  return(M)
}

B_index_from_Mraw_index = function(M,R,N) {
 B = (M - (N-1)/R )/N
  return(B)
}
