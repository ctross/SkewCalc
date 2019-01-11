#' Calculate the M index from the B index and relevant data.
#'
#' @param B The B index value estimated from a dataset.
#' @param R Total number of RS events in the sample from which B was estimated.
#' @param R Total number of individuals in the sample from which B was estimated.
#' @return The M index.
#' @examples
#' B = 0.1 
#' R = 50
#' N = 10
#' M_index_from_B_index(B, R, N)

M_index_from_B_index = function(B, R, N) {
  M = sqrt(B*R - 1/N + 1)
  return(M)
}
