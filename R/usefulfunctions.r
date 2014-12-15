##################################################################################
flist_untag  <- function (flist)
{
    for (i in 1:length(flist)) {
        if (class(flist[[i]]) == "<-") {
            flist[[i]][[1]] <- as.name("~")
        }
        flist[[i]] <- eval(flist[[i]])
    }
    as.list(flist)
}

################################################################################# Define Link Function
logistic<-function (yy)
{
    p <- 1/(1 + exp(-yy))
    p <- ifelse(yy == Inf, 1, p)
    p
}

################################################################################### Convergence Check
xcheckconvergence <-function (model)
{
    the.class <- class(model)[1]
    k <- 0
    if (the.class == "mle2") {
        if (model@details$convergence != 0) {
            k <- model@details$convergence
        }
    }
    if (the.class == "map") {
        if (model@optim$convergence != 0) {
            k <- model@optim$convergence
        }
    }
    if (k > 0) {
        message(paste("Caution, model may not have converged."))
        if (k == 1) {
            message("Code 1: Maximum iterations reached.")
        }
        if (k == 10) {
            message("Code 10: Degenerate Nelder-Mead simplex.")
        }
    }
}
###################################################################### Confidence Intervals
PCI<-function (samples, prob = P)
{
    a <- (1 - prob)/2
    quantile(samples, probs = c(a, 1 - a))
}

HPDI<-function (samples, prob = P)
{
    class.samples <- class(samples)[1]
    coerce.list <- c("numeric", "matrix", "data.frame", "integer",
        "array")
    if (class.samples %in% coerce.list) {
        samples <- as.mcmc(samples)
    }
    x <- coda::HPDinterval(samples, prob = prob)
    result <- c(x[1], x[2])
    names(result) <- c(paste("lower", prob), paste("upper", prob))
    result
}


 ######################################### Calc Skew - 
M_Mc<-function(Pred) {
NN<-length(Pred)/2
Si<-c()
M<-c()

	K <- sum(Pred[1:NN]);
	F <- sum(Pred[(NN+1):(2*NN)]);
	Si <- ((Pred[1:NN]/K)-(Pred[(NN+1):(2*NN)]/F)) * ((Pred[1:NN]/K)-(Pred[(NN+1):(2*NN)]/F)) ;
	S <- sum(Si);
 	M[1] <- sqrt(NN * S);
	M[2] <- sqrt(NN * S)*sqrt(mean(Pred[1:NN]));
  return(M)
}

 ######################################### Calc Skew 
M_Mc<-function(Pred) {
NN<-length(Pred)/2
Si<-c()
M<-c()

	K <- sum(Pred[1:NN]);
	F <- sum(Pred[(NN+1):(2*NN)]);
	Si <- ((Pred[1:NN]/K)-(Pred[(NN+1):(2*NN)]/F)) * ((Pred[1:NN]/K)-(Pred[(NN+1):(2*NN)]/F)) ;
	S <- sum(Si);
 	M[1] <- sqrt(NN * S);
	M[2] <- sqrt(NN * S)*sqrt(mean(Pred[1:NN]));
  return(M)
}

########################################################### Sample Indices

B_Index <- function(ki,ni) {
	Nt <- sum(ni)
	Nbar <- Nt / max(ni)
	K <- sum(ki)
	if(K>0){
		B = sum((ki / K - ni / Nt)^2) - (1 / K) * (1 - 1 / Nbar)
	} else {
		B = 0
	}
	B
}

M_Index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	C <- sqrt(N * S)
	C
}

M_Index_From_B_Index <- function(N,B,K) {
	C <- sqrt(N * B + (N-1)/K)
	C
}

Gini<-function (x, corr = FALSE, na.rm = TRUE) 
{
    if (!na.rm && any(is.na(x))) 
        return(NA_real_)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1L:n)
    G <- 2 * G/sum(x) - (n + 1L)
    if (corr) 
        G/(n - 1L)
    else G/n
}

Mc_Index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	D1 <- sqrt(N * S)*sqrt(mean(ki))
  D1
}




