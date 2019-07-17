

function (samples, prob = 0.89) 
{
    x <- sapply(prob, function(p) {
        a <- (1 - p)/2
        quantile(samples, probs = c(a, 1 - a))
    })
    n <- length(prob)
    result <- rep(0, n * 2)
    for (i in 1:n) {
        low_idx <- n + 1 - i
        up_idx <- n + i
        result[low_idx] <- x[1, i]
        result[up_idx] <- x[2, i]
        a <- (1 - prob[i])/2
        names(result)[low_idx] <- concat(round(a * 100, 0), "%")
        names(result)[up_idx] <- concat(round((1 - a) * 100, 
            0), "%")
    }
    return(result)
}
