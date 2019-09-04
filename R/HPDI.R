#' HPDI - Stolen from rethinking by R. McElreath
#'
#' @param samples A vector of data-points.
#' @param prob A percentile interval.
#' @return The higest posterior density interval.
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' HPDI(RS, 0.9)

HPDI <- function (samples, prob = 0.89) 
{
    class.samples <- class(samples)[1]
    coerce.list <- c("numeric", "matrix", "data.frame", "integer", 
        "array")
    if (class.samples %in% coerce.list) {
        samples <- coda::as.mcmc(samples)
    }
    x <- sapply(prob, function(p) coda::HPDinterval(samples, 
        prob = p))
    n <- length(prob)
    result <- rep(0, n * 2)
    for (i in 1:n) {
        low_idx <- n + 1 - i
        up_idx <- n + i
        result[low_idx] <- x[1, i]
        result[up_idx] <- x[2, i]
        names(result)[low_idx] <- concat("|", prob[i])
        names(result)[up_idx] <- concat(prob[i], "|")
    }
    return(result)
}
