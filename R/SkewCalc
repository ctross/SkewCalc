
# Code and Models By Cody T. Ross, SFI, UC Davis
# Many functions used are taken from the Rethinking Package from Richard McElreath

################################################################################### Read in Maximum A Posteriori Estimating Algorithm
map<-function (flist, data, start, method = "BFGS", hessian = TRUE,
    debug = FALSE, verbose = FALSE, ...)
{
##################### Required Dependencies
    library(MASS)
    library(bbmle)
    library(coda)
# map function is taken directly from Richard McElreath, Rethinking package
    if (missing(flist))
        stop("Formula required.")
    if (class(flist) != "list") {
        if (class(flist) == "formula") {
            flist <- list(flist)
        }
        else {
            stop("Formula or list of formulas required.")
        }
    }
    if (missing(start))
        start <- list()
    if (missing(data))
        stop("'data' required.")
    if (!(class(data) %in% c("list", "data.frame"))) {
        stop("'data' must be of class list or data.frame.")
    }
    flist.orig <- flist
    flist <- flist_untag(flist)
    density_sub_list <- list(normal = "dnorm", binomial = "dbinom",
        poisson = "dpois", gamma = "dgamma2", betabinomial = "dbetabinom",
        gammapoisson = "dgampois", cauchy = "dcauchy", uniform = "dunif",
        laplace = "dlaplace")
    idx_marker_string <- "_._"
    sample_from_prior <- function(f) {
        RHS <- f[[3]]
        the_density <- as.character(RHS[[1]])
        the_rdensity <- the_density
        substr(the_rdensity, 1, 1) <- "r"
        pars <- vector(mode = "list", length = length(RHS))
        pars[[1]] <- 1
        for (i in 1:(length(pars) - 1)) {
            pars[[i + 1]] <- RHS[[i + 1]]
            pars[[i + 1]] <- eval(pars[[i + 1]], as.list(data))
        }
        result <- do.call(the_rdensity, args = pars)
        return(result)
    }
    dparser <- function(flist, e) {
        r <- sapply(flist, function(i) sum(eval(parse(text = i),
            envir = e)))
        -sum(r)
    }
    pars_to_vectors <- function(pars_orig, veclist) {
        pars_new <- list()
        for (i in 1:length(pars_orig)) {
            x <- strsplit(names(pars_orig)[i], idx_marker_string)[[1]]
            if (length(x) == 1) {
                pars_new[[x]] <- pars_orig[[i]]
            }
        }
        for (i in 1:length(veclist)) {
            newvec <- rep(0, veclist[[i]]$n)
            for (j in 1:veclist[[i]]$n) {
                name_orig <- paste(veclist[[i]]$name, idx_marker_string,
                  j, sep = "", concat = "")
                newvec[j] <- pars_orig[[name_orig]]
            }
            pars_new[[veclist[[i]]$name]] <- newvec
        }
        return(pars_new)
    }
    make_minuslogl <- function(pars, flist, data, veclist) {
        if (length(veclist) > 0)
            pars <- pars_to_vectors(pars, veclist)
        e <- list(as.list(data), as.list(pars))
        e <- unlist(e, recursive = FALSE)
        dparser(flist, e)
    }
    link.names <- c("log", "logit")
    invlink.names <- c("exp", "logistic")
    formula2text <- function(f) {
        RHS <- f[[3]]
        LHS <- f[[2]]
        fname <- as.character(RHS[[1]])
        if (fname == "+" | fname == "*" | fname == "-" | fname ==
            "/" | fname == "%*%" | fname %in% invlink.names) {
            thetext <- list(as.character(LHS), paste(deparse(RHS),
                collapse = " "))
        }
        else {
            n_args <- length(RHS)
            args_list <- as.list(RHS)
            if (class(LHS) == "call") {
                if (as.character(LHS[[1]]) == "[") {
                  ival <- suppressWarnings(as.numeric(as.character(LHS[[3]])))
                  if (is.na(ival)) {
                    LHS <- as.character(LHS[[2]])
                  }
                  else {
                    LHS <- deparse(LHS)
                  }
                }
            }
            args_list[[1]] <- LHS
            args_list[[n_args + 1]] <- "TRUE"
            args_names <- names(RHS)[-1]
            if (is.null(args_names))
                args_names <- rep("", n_args - 1)
            args_names <- c("x", args_names, "log")
            for (i in 1:length(args_names)) {
                if (args_names[i] == "") {
                  args_names[i] <- args_list[i]
                }
                else {
                  args_names[i] <- paste(args_names[i], args_list[i],
                    sep = "=")
                }
            }
            args_text <- paste(args_names, collapse = " , ")
            thetext <- paste(RHS[[1]], "(", args_text, ")", sep = "")
        }
        return(thetext)
    }
    mygrep <- function(target, replacement, x, add.par = TRUE) {
        wild <- "[()=*+ ]"
        pattern <- paste(wild, target, wild, sep = "", collapse = "")
        m <- regexpr(pattern, x)
        if (m == -1)
            return(x)
        s <- regmatches(x = x, m = m)
        if (add.par == TRUE)
            replacement <- paste("(", replacement, ")", collapse = "")
        w.start <- substr(s, 1, 1)
        w.end <- substr(s, nchar(s), nchar(s))
        r <- paste(w.start, replacement, w.end, sep = "", collapse = "")
        gsub(pattern = s, replacement = r, x = x, fixed = TRUE)
    }
    pars_with_priors <- list()
    if (length(flist) > 1) {
        flag_flatten <- FALSE
        for (i in 2:length(flist)) {
            if (!(class(flist[[i]]) == "formula"))
                stop("Input not a formula.")
            LHS <- flist[[i]][[2]]
            if (class(LHS) == "call") {
                fname <- as.character(LHS[[1]])
                if (fname == "c" | fname == "[" | fname %in%
                  link.names) {
                  if (fname == "c") {
                    newflist <- list()
                    num_pars <- length(LHS) - 1
                    for (j in 1:num_pars) {
                      newflist[[j]] <- flist[[i]]
                      newflist[[j]][[2]] <- LHS[[j + 1]]
                      pars_with_priors[[as.character(LHS[[j +
                        1]])]] <- 1
                    }
                    flist[[i]] <- newflist
                    flag_flatten <- TRUE
                  }
                  if (fname %in% link.names) {
                    the.link <- as.character(LHS[[1]])
                    flist[[i]][[2]] <- flist[[i]][[2]][[2]]
                    the.invlink <- invlink.names[which(link.names ==
                      the.link)]
                    old.RHS <- flist[[i]][[3]]
                    flist[[i]][[3]] <- as.call(list(as.name(the.invlink),
                      old.RHS))
                  }
                  if (fname == "[") {
                    pars_with_priors[[deparse(LHS)]] <- 1
                  }
                }
                else {
                  stop(paste("Invalid prior specification:",
                    deparse(flist[[i]])))
                }
            }
            else {
                pars_with_priors[[as.character(LHS)]] <- 1
            }
        }
        if (flag_flatten)
            flist <- unlist(flist, recursive = FALSE)
    }
    if (debug)
        print(flist)
    flist2 <- lapply(flist, formula2text)
    if (debug)
        print(flist2)
    links <- list()
    if (length(flist2) > 1) {
        for (i in length(flist2):2) {
            if (class(flist2[[i]]) == "list") {
                LHS <- flist2[[i]][[1]]
                RHS <- flist2[[i]][[2]]
                lik_save <- flist2[[1]]
                flist2[[1]] <- mygrep(LHS, RHS, flist2[[1]],
                  add.par = FALSE)
                if (flist2[[1]] != lik_save) {
                  links[[length(links) + 1]] <- flist2[[i]]
                }
                if (i > 2) {
                  for (j in (i - 1):2) {
                    if (class(flist2[[j]]) == "list") {
                      flist2[[j]][[2]] <- mygrep(LHS, RHS, flist2[[j]][[2]],
                        add.par = TRUE)
                    }
                  }
                }
                pars_with_priors[[LHS]] <- NULL
            }
        }
        flist3 <- list()
        j <- 1
        for (i in 1:length(flist2)) {
            if (class(flist2[[i]]) != "list") {
                flist3[[j]] <- flist2[[i]]
                j <- j + 1
            }
        }
        flist2 <- flist3
    }
    flist.ll <- flist2[[1]]
    if (debug) {
        print(flist.ll)
    }
    if (debug) {
        print(start)
        print(pars_with_priors)
    }
    pars <- start
    if (any(!(names(pars_with_priors) %in% names(pars)))) {
        bad_pars <- names(pars_with_priors)[!(names(pars_with_priors) %in%
            names(pars))]
        if (verbose == TRUE)
            message(paste("Sampling start values from priors for:",
                paste(bad_pars, collapse = " ")))
        for (k in bad_pars) {
            for (g in 2:length(flist)) {
                if (class(flist[[g]][[2]]) == "call") {
                  the_par_with_index <- deparse(flist[[g]][[2]])
                  the_par <- as.character(flist[[g]][[2]][[2]])
                  the_index_or_var <- as.character(flist[[g]][[2]][[3]])
                  if (the_par_with_index == k) {
                    idx_num <- suppressWarnings(as.numeric(the_index_or_var))
                    if (!is.na(idx_num)) {
                      if (is.null(start[[the_par]])) {
                        max_index <- 1
                        for (h in 2:length(flist)) {
                          if (class(flist[[h]][[2]]) == "call") {
                            if (as.character(flist[[h]][[2]][[2]]) ==
                              "[" & as.character(flist[[h]][[2]][[3]]) ==
                              the_par) {
                              nval <- suppressWarnings(as.numeric(flist[[h]][[2]][[3]]))
                              if (!is.null(nval)) {
                                max_index <- max(max_index, nval)
                              }
                            }
                          }
                        }
                        start[[the_par]] <- rep(0, max_index)
                        the_index <- as.numeric(flist[[g]][[2]][[3]])
                        start[[the_par]][the_index] <- sample_from_prior(flist[[g]])
                      }
                      else {
                        the_index <- as.numeric(flist[[g]][[2]][[3]])
                        start[[the_par]][the_index] <- sample_from_prior(flist[[g]])
                      }
                    }
                    if (is.na(idx_num)) {
                      the_var <- as.character(flist[[g]][[2]][[3]])
                      n_unique <- length(unique(data[[the_var]]))
                      start[[the_par]] <- replicate(n_unique,
                        sample_from_prior(flist[[g]]))
                    }
                  }
                }
                else {
                  the_par <- paste(as.character(flist[[g]][[2]]),
                    collapse = "")
                  if (the_par == k) {
                    f <- flist[[g]]
                    theta <- sample_from_prior(f)
                    start[[k]] <- theta
                  }
                }
            }
        }
    }
    if (FALSE) {
        if (any(!(names(pars) %in% names(pars_with_priors)))) {
            flat_pars <- names(pars)[!(names(pars) %in% names(pars_with_priors))]
            message(paste("Using flat priors for:", paste(flat_pars,
                collapse = " ")))
        }
    }
    pars <- start
    pars_flat <- list()
    veclist <- list()
    for (i in 1:length(pars)) {
        n <- length(pars[[i]])
        par_name <- names(pars)[i]
        if (n == 1) {
            pars_flat[[par_name]] <- pars[[i]]
        }
        else {
            for (j in 1:n) {
                new_name <- concat(par_name, idx_marker_string,
                  j)
                pars_flat[[new_name]] <- pars[[i]][j]
            }
            veclist[[par_name]] <- list(name = par_name, n = n)
        }
    }
    pars <- pars_flat
    fit <- try(suppressWarnings(optim(par = pars, fn = make_minuslogl,
        flist = flist2, data = data, veclist = veclist, hessian = hessian,
        method = method, ...)), silent = TRUE)
    if (class(fit) == "try-error") {
        msg <- attr(fit, "condition")$message
        objnotfound <- grep("object '.*' not found", msg)
        if (length(objnotfound) > 0) {
            obj_name <- regmatches(msg, gregexpr("'.*'", msg))
            out_msg <- paste("Cannot find ", obj_name, ".\nIf this is a parameter, try defining a prior for it or providing a start value.\nIf this is a variable, make sure it is in the data list.",
                collapse = "", sep = "")
            stop(out_msg)
        }
        objnotfound <- grep("initial value in 'vmmin' is not finite",
            msg)
        if (length(objnotfound) > 0) {
            out_msg <- paste(msg, "\nThe start values for the parameters were invalid. This could be caused by missing values (NA) in the data or by start values outside the parameter constraints. If there are no NA values in the data, try using explicit start values.",
                collapse = "", sep = "")
            stop(out_msg)
        }
        objnotfound <- grep("non-finite finite-difference value",
            msg)
        if (length(objnotfound) > 0) {
            out_msg <- paste(msg, "\nStart values for parameters may be too far from MAP.\nTry better priors or use explicit start values.\nIf you sampled random start values, just trying again may work.\nStart values used in this attempt:\n",
                paste(names(start), "=", start, collapse = "\n"),
                collapse = "", sep = "")
            stop(out_msg)
        }
        stop(msg)
    }
    if (hessian) {
        vcov <- try(solve(fit$hessian))
        if (class(vcov) == "try-error") {
            warning("Error when computing variance-covariance matrix (Hessian). Fit may not be reliable.")
            vcov <- matrix(NA, nrow = length(pars), ncol = length(pars))
        }
    }
    else {
        vcov <- matrix(NA, nrow = length(pars), ncol = length(pars))
    }
    fit$minuslogl <- make_minuslogl(fit$par, flist = flist.ll,
        data = data, veclist = veclist)
    fmll <- function(pars) make_minuslogl(pars, flist = flist.ll,
        data = data, veclist = veclist)
    coefs <- fit$par
    for (i in 1:length(coefs)) {
        a_split <- strsplit(names(coefs)[i], idx_marker_string)[[1]]
        if (length(a_split) > 1) {
            new_name <- concat(a_split[1], "[", a_split[2], "]")
            names(coefs)[i] <- new_name
        }
    }
    m <- new("map", call = match.call(), coef = coefs, vcov = vcov,
        optim = fit, data = as.list(data), start = start, formula = flist.orig,
        formula_parsed = flist2, fminuslogl = fmll, links = links)
    attr(m, "df") <- length(m@coef)
    attr(m, "veclist") <- veclist
    if (!missing(data))
        attr(m, "nobs") = length(data[[1]])
    xcheckconvergence(m)
    m
}

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


################################################################################### Create Class
setClass("map", representation( call = "language",
coef = "numeric",
vcov = "matrix",
optim = "list",
data = "list",
start = "list",
formula = "list",
formula_parsed = "list",
fminuslogl = "function",
links = "list"))
setMethod("coef", "map", function(object) {
object@coef
})
setMethod("vcov", "map", function (object, ...) { object@vcov } )
setMethod("nobs", "map", function (object, ...) { attr(object,"nobs") } )
setGeneric("stancode",
function( object , ... ) {
print(class(object))
}
)
setMethod("stancode", "map",
function(object) {
# compile Stan code through map2stan
m <- map2stan( object@formula , data=object@data , start=object@start , sample=FALSE )
# display
cat( m$model )
return( invisible( m$model ) )
}
)
setMethod("logLik", "map",
function (object, ...)
{
if(length(list(...)))
warning("extra arguments discarded")
val <- -object@optim[['minuslogl']]
attr(val, "df") <- length(object@coef)
attr(val, "nobs") <- attr(object,"nobs")
class(val) <- "logLik"
val
})
setMethod("deviance", "map",
function (object, ...)
{
-2*logLik(object)
})
setMethod("AIC", "map",
function (object, ...)
{
-2*logLik(object) + 2*length(coef(object))
})
setMethod("show", "map", function(object){
cat("\nMaximum a posteriori (MAP) model fit\n")
cat("\nFormula:\n")
for ( i in 1:length(object@formula) ) {
print( object@formula[[i]] )
}
cat("\nMAP values:\n")
print(coef(object))
cat("\nLog-likelihood: ")
cat(round(as.numeric(logLik(object)),2),"\n")
if ( object@optim$convergence > 0 )
cat("\nWarning: optimization did not converge (code ",
object@optim$convergence, ": " , object@optim$message,")\n",sep="")
})
setMethod("summary", "map", function(object){
precis(object)
})
setGeneric("extract.samples",
function( object , n=10000 , clean.names=TRUE , ... ) {
require(MASS)
mu <- 0
if ( class(object)[1] %in% c("mer","bmer","glmerMod","lmerMod") ) {
mu <- fixef(object)
} else {
mu <- xcoef(object)
}
result <- as.data.frame( mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
if ( clean.names==TRUE ) {
# convert (Intercept) to Intercept
for ( i in 1:ncol(result) ) {
if ( colnames(result)[i] == "(Intercept)" ) {
colnames(result)[i] <- "Intercept"
}
}
}
result
}
)
setMethod("extract.samples", "map",
function(object,n=1e4,...){
require(MASS)
mu <- object@coef
result <- as.data.frame( mvrnorm( n=n , mu=mu , Sigma=vcov(object) ) )
# convert vector parameters to vectors in list
veclist <- attr(object,"veclist")
name_head <- function(aname) strsplit( aname , "[" , fixed=TRUE )[[1]][1]
name_index <- function(aname) as.numeric(regmatches( aname , regexec( "\\[(.+)\\]" , aname ) )[[1]][2])
if ( length(veclist) > 0 ) {
new_result <- list()
# copy non-vector samples into new list
for ( i in 1:length(result) ) {
if ( !( name_head(names(result)[i]) %in% names(veclist) ) ) {
new_result[[ names(result)[i] ]] <- result[[i]]
}
}#i
# now build vectors out of paramters with [n] in name
for ( i in 1:length(veclist) ) {
# empty matrix with parameters on cols and samples on rows
# so n-by-m, where m in number of pars in vector and n is number of samples
new_matrix <- matrix( 0 , ncol=veclist[[i]]$n , nrow=n )
for ( j in 1:length(result) ) {
if ( name_head(names(result)[j]) == names(veclist)[i] ) {
the_index <- name_index( names(result)[j] )
new_matrix[,the_index] <- result[[j]]
}
}
new_result[[ names(veclist)[i] ]] <- new_matrix
}#i
result <- new_result
}
# return result
result
})

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
###################################################################### Intervals
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
#################################################################################### Precis
precis.whitelist <- data.frame(
class=c("map","map2stan","lm","glm","mle2","mer","bmer","polr","data.frame","clmm","clmm2","list","stanfit","lmerMod","glmerMod") ,
coef.method=c("coef","coef","coef","coef","coef","fixef.plus","fixef.plus","polr","chain","coef","coef","mcarray","stanfit","fixef.plus","fixef.plus") ,
vcov.method=c("vcov","vcov","vcov","vcov","vcov","vcov.VarCorr","vcov.VarCorr","vcov","chain","vcov","vcov","mcarray","stanfit","vcov.VarCorr","vcov.VarCorr") ,
nobs.method=c("nobs","nobs","nobs","nobs","mle2","mer","mer","nobs","chain","nobs","nobs","chain","stanfit","mer","mer")
)
# precis class definition and show method
setClass( "precis" , representation( output="data.frame" , digits="numeric" ) )
precis.show <- function( object ) {
#print( round( object@output , object@digits ) )
r <- format_show( object@output , digits=c('default__'=object@digits,'n_eff'=0) )
print(r)
}
setMethod( "show" , "precis" , function(object) precis.show(object) )
precis.plot <- function( x , y , pars , col.ci="black" , xlab="Value" , ... ) {
x <- x@output
if ( !missing(pars) ) {
x <- x[pars,]
}
n <- nrow(x)
mu <- x[n:1,1]
left <- x[[3]][n:1]
right <- x[[4]][n:1]
set_nice_margins()
dotchart( mu , labels=rownames(x)[n:1] , xlab=xlab , xlim=c(min(left),max(right)) , ... )
for ( i in 1:length(mu) ) lines( c(left[i],right[i]) , c(i,i) , lwd=2 , col=col.ci )
abline( v=0 , lty=1 , col=col.alpha("black",0.15) )
}
setMethod( "plot" , "precis" , function(x,y,...) precis.plot(x,y,...) )
# function to process a list of posterior samples from extract.samples into a summary table
# needed because as.data.frame borks the ordering of matrix parameters like varying effects
postlistprecis <- function( post , prob=0.95 ) {
n_pars <- length(post)
result <- data.frame( Mean=0 , StdDev=0 , lower=0 , upper=0 )
r <- 1
for ( k in 1:n_pars ) {
dims <- dim( post[[k]] )
if ( length(dims)==1 ) {
# single parameter
hpd <- as.numeric( HPDI( post[[k]] , prob=prob ) )
result[r,] <- c( mean(post[[k]]) , sd(post[[k]]) , hpd[1] , hpd[2] )
rownames(result)[r] <- names(post)[k]
r <- r + 1
}
if ( length(dims)==2 ) {
# vector of parameters
# loop over
for ( i in 1:dims[2] ) {
hpd <- as.numeric( HPDI( post[[k]][,i] , prob=prob ) )
result[r,] <- c( mean(post[[k]][,i]) , sd(post[[k]][,i]) , hpd[1] , hpd[2] )
rownames(result)[r] <- concat( names(post)[k] , "[" , i , "]" )
r <- r + 1
}
}
if ( length(dims)==3 ) {
# matrix of parameters
for ( i in 1:dims[2] ) {
for ( j in 1:dims[3] ) {
hpd <- as.numeric( HPDI( post[[k]][,i,j] , prob=prob ) )
result[r,] <- c( mean(post[[k]][,i,j]) , sd(post[[k]][,i,j]) , hpd[1] , hpd[2] )
rownames(result)[r] <- concat( names(post)[k] , "[" , i , "," , j , "]" )
r <- r + 1
}
}
}
}
colnames(result)[3:4] <- c(paste("lower", prob), paste("upper", prob))
result
}
precis <- function( model , depth=1 , pars , ci=TRUE , level=0.95 , corr=FALSE , digits=2 , warn=TRUE ) {
the.class <- class(model)[1]
found.class <- FALSE
if ( the.class=="numeric" ) {
# single vector of values
# coerce to data frame
model <- as.data.frame(model)
the.class <- class(model)[1]
}
if ( any( precis.whitelist$class==the.class ) ) found.class <- TRUE
if ( the.class=="list" )
if ( class( model[[1]] ) != "mcarray" ) found.class <- FALSE
if ( found.class==TRUE ) {
est <- xcoef( model )
se <- xse( model )
if ( corr==TRUE ) Rho <- xrho( model )
}
if ( found.class==FALSE ) {
message( paste("No handler found for model of class",the.class) )
return(invisible())
}
# format
result <- data.frame( est=est , se=se )
colnames(result) <- c("Mean","StdDev")
if ( ci==TRUE ) {
ci <- confint.quad( est=est , se=se , level=level )
if ( the.class=="data.frame" ) {
# HPDI from samples
ci <- t( apply( model , 2 , HPDI , prob=level ) )
}
result <- cbind( result , ci )
if ( the.class=="map2stan" ) {
# HPDI from samples
post <- extract.samples(model)
result <- postlistprecis( post , prob=level )
}
if ( the.class=="stanfit" ) {
# HPDI from samples
post <- extract.samples(model)
post[['lp__']] <- NULL
result <- postlistprecis( post , prob=level )
}
}
if ( the.class=="map2stan" | the.class=="stanfit" ) {
# add n_eff to result
require(rstan)
if ( the.class=="map2stan" )
the_summary <- summary( model@stanfit )$summary
else
the_summary <- summary( model )$summary
n_eff <- the_summary[,'n_eff']
n_eff <- n_eff[ -which(names(n_eff)=="lp__") ]
Rhat <- the_summary[,'Rhat']
Rhat <- Rhat[ -which(names(Rhat)=="lp__") ]
if ( the.class=="map2stan" ) {
n_eff <- n_eff[ -which(names(n_eff)=="dev") ]
Rhat <- Rhat[ -which(names(Rhat)=="dev") ]
}
result <- cbind( result , n_eff , Rhat )
}
if ( corr==TRUE ) {
result <- cbind( result , Rho )
}
#if ( type.s==TRUE )
# result[,"Pr(S)"] <- format.pval( type.s( est , se ) )
if ( precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]=="vcov.VarCorr" ) {
message( "Quadratic approximation (standard errors) unreliable for variance components. Use MCMC to estimate precision of variance components." )
}
# deal with depth
if ( depth==1 ) {
hits <- regexpr("]",rownames(result),fixed=TRUE)
hits_idx <- which( hits > -1 )
if ( length(hits_idx)>0 ) {
result <- result[-hits_idx,]
message( paste( length(hits_idx) , "vector or matrix parameters omitted in display. Use depth=2 to show them." ) )
}
}
# deal with pars list
if ( !missing(pars) ) {
# have to handle vector/matrix parameters,
# so split off any [.] and then prune with names only
clean_names <- as.character( sapply( rownames(result) , function(s) strsplit( s , "[" , fixed=TRUE )[[1]][1] ) )
result <- result[ clean_names %in% pars , ]
}
# result
new( "precis" , output=result , digits=digits )
}
####
xcoef <- function( model ) {
the.class <- class(model)[1]
the.method <- precis.whitelist$coef.method[ precis.whitelist$class==the.class ]
if ( the.method=="coef" ) {
result <- coef(model)
}
if ( the.method=="fixef" ) {
result <- fixef(model)
}
if ( the.method=="polr" ) {
result <- summary(model)$coefficients[,1]
}
if ( the.method=="chain" ) {
# average of chains
result <- apply( model , 2 , mean )
}
if ( the.method=="stanfit" ) {
result <- summary( model )$summary[,1]
}
if ( the.method=="mcarray" ) {
# jags.samples result, hopefully
result <- NULL
result.names <- NULL
for ( j in 1:length(model) ) {
# explode compact arrays of coefficients
dims <- dim( model[[j]] )
if ( length(dims)==3 ) {
est <- rep(0,dims[1])
for ( k in 1:dims[1] ) {
est[k] <- mean( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
}
result <- c( result , est )
if ( dims[1] > 1 ) {
newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
} else {
newnames <- names(model)[j]
}
result.names <- c( result.names , newnames )
} # if dim 3
} # for each array
names(result) <- result.names
}
if ( the.method=="fixef.plus" ) {
# fixef from lmer plus variance components
result <- fixef(model)
vc <- VarCorr(model)
clusters <- names(vc)
for( i in 1:length(clusters) ) {
sigma <- sqrt(diag(vc[[i]]))
names(sigma) <- paste( rownames(vc[[i]]) , clusters[i] , sep="|" )
names(sigma) <- paste( "(" , names(sigma) , ")" , sep="" )
result <- c( result , sigma )
}
sigma.resid <- attr( vc , "sc" )
if ( !is.na(sigma.resid) ) result['(residual)'] <- sigma.resid
}
xcheckconvergence( model )
result
}
xcheckconvergence <- function( model ) {
the.class <- class(model)[1]
k <- 0
if ( the.class=="mle2" ) {
if ( model@details$convergence != 0 ) {
k <- model@details$convergence
}
}
if ( the.class=="map" ) {
if ( model@optim$convergence != 0 ) {
k <- model@optim$convergence
}
}
if ( k > 0 ) {
message( paste("Caution, model may not have converged.") )
if ( k==1 ) {
message( "Code 1: Maximum iterations reached." )
}
if ( k==10 ) {
message( "Code 10: Degenerate Nelder-Mead simplex." )
}
}
}
xse <- function( model ) {
the.class <- class(model)[1]
the.method <- precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]
if ( the.method=="vcov" ) {
result <- sqrt(diag(vcov(model)))
}
if ( the.method=="vcov.VarCorr" ) {
result <- sqrt(diag( as.matrix(vcov(model)) ))
num.sigma <- length( xcoef(model) ) - length( fixef(model) )
result <- c( result , rep(NA,num.sigma) )
}
if ( the.method=="chain" ) {
# sd of chains
result <- apply( model , 2 , sd )
}
if ( the.method=="stanfit" ) {
result <- summary( model )$summary[,3]
}
if ( the.method=="mcarray" ) {
# jags.samples result, hopefully
result <- NULL
result.names <- NULL
for ( j in 1:length(model) ) {
# explode compact arrays of coefficients
dims <- dim( model[[j]] )
if ( length(dims)==3 ) {
est <- rep(0,dims[1])
for ( k in 1:dims[1] ) {
est[k] <- sd( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
}
result <- c( result , est )
if ( dims[1] > 1 ) {
newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
} else {
newnames <- names(model)[j]
}
result.names <- c( result.names , newnames )
} # if dim 3
} # for each array
names(result) <- result.names
}
result
}
xrho <- function( model ) {
the.class <- class(model)[1]
the.method <- precis.whitelist$vcov.method[ precis.whitelist$class==the.class ]
if ( the.method=="vcov" ) {
result <- cov2cor( vcov(model) )
}
if ( the.method=="vcov.VarCorr" ) {
result <- sqrt(diag( as.matrix(vcov(model)) ))
num.sigma <- length( xcoef(model) ) - length( fixef(model) )
result <- c( result , rep(NA,num.sigma) )
}
if ( the.method=="chain" ) {
# sd of chains
result <- apply( model , 2 , sd )
}
if ( the.method=="stanfit" ) {
result <- summary( model )$summary[,3]
}
if ( the.method=="mcarray" ) {
# jags.samples result, hopefully
result <- NULL
result.names <- NULL
for ( j in 1:length(model) ) {
# explode compact arrays of coefficients
dims <- dim( model[[j]] )
if ( length(dims)==3 ) {
est <- rep(0,dims[1])
for ( k in 1:dims[1] ) {
est[k] <- sd( as.vector(model[[j]][ k , , ]) ) # marginalize over iteration and chain
}
result <- c( result , est )
if ( dims[1] > 1 ) {
newnames <- paste( names(model)[j] , "[" , 1:dims[1] , "]" , sep="" )
} else {
newnames <- names(model)[j]
}
result.names <- c( result.names , newnames )
} # if dim 3
} # for each array
names(result) <- result.names
}
result
}
xnobs <- function( model ) {
the.class <- class(model)[1]
the.method <- precis.whitelist$nobs.method[ precis.whitelist$class==the.class ]
if ( the.method=="nobs" ) result <- nobs(model)
if ( the.method=="mle2" ) result <- length(model@data[[1]])
if ( the.method=="mer" ) result <- nrow(model@frame)
if ( the.method=="chain" ) result <- nrow(model)
if ( the.method=="stanfit" ) result <- 0
result
}
# row-by-row matrix formatting function
rrformat <- function( matrix , digits=2 , width=7 ) {
if ( length(digits)==1 ) digits <- rep(digits,nrow(matrix))
result <- matrix
for ( i in 1:nrow(matrix) ) {
result[i,] <- format( round(matrix[i,],digits[i]) , width=width )
}
result
}

confint.quad<-function (model = NULL, est, se, level = 0.95)
{
    if (!is.null(model)) {
        found.class <- FALSE
        if (class(model) == "lm") {
            est <- coef(model)
            se <- summary(model)$coef[, 2]
            found.class <- TRUE
        }
        if (class(model) == "mle2") {
            est <- coef(model)
            se <- summary(model)@coef[, 2]
            found.class <- TRUE
        }
        if (found.class == FALSE) {
            return(paste("Cannot find handler for model of class",
                class(model)))
        }
    }
    n <- length(est)
    mat <- matrix(c(rep(-1, n), rep(1, n)), nrow = n)
    p <- (1 - level)/2
    z <- -qnorm(p)
    ci <- est + mat * (se * z)
    rownames(ci) <- names(est)
    lowlab <- paste(format(p * 100, nsmall = 1), "%", sep = "")
    hilab <- paste(format((1 - p) * 100, nsmall = 1), "%", sep = "")
    colnames(ci) <- c(lowlab, hilab)
    ci
}

 format_show <- function( x , digits ) {
r <- as.data.frame(lapply( 1:length(x) ,
function(i) {
if ( names(x)[i] %in% names(digits) )
round( x[[i]] , digits[names(x)[i]] )
else
round(x[[i]], digits['default__'] );
} ) )
names(r) <- names(x)
rownames(r) <- rownames(x)
return(r)
}


############## Load Function to Process Posterior
 sample.qa.posterior <- function (model, n = Samples, clean.names = TRUE, model.weights = "AICc",
    nobs = 0, add.names = FALSE, fill.na = 0, verbose = FALSE){
    require(MASS)
    require(bbmle)
    getdf <- function(x) {
        if (!is.null(df <- attr(x, "df")))
            return(df)
        else if (!is.null(df <- attr(logLik(x), "df")))
            return(df)
    }
    myBIC <- function(x, nobs) {
        k <- getdf(x)
        as.numeric(-2 * logLik(x) + log(nobs) * k)
    }
    if (class(model)[1] == "list") {
        if (length(model) < 2) {
            return(sample.qa.posterior(model[[1]], n = n, nobs = nobs))
        }
        valid.methods <- c("AIC", "AICc", "BIC")
        use.custom <- FALSE
        if (class(model.weights) == "numeric") {
            if (length(model.weights) != length(model)) {
                stop("Custom model weights must be same length as list of models.")
            }
            else {
                use.custom <- TRUE
                post <- model.weights
            }
        }
        else {
            if (!any(model.weights == valid.methods)) {
                stop(paste("Unknown model averaging method:",
                  model.weights))
            }
        }
        if (use.custom == FALSE) {
            if (model.weights == "AIC")
                factors <- sapply(model, AIC)
            if (nobs == 0) {
                if (model.weights == "AICc")
                  factors <- sapply(model, AICc)
                if (model.weights == "BIC")
                  factors <- sapply(model, BIC)
            }
            else {
                if (model.weights == "AICc")
                  factors <- sapply(model, function(z) AICc(z,
                    nobs = nobs))
                if (model.weights == "BIC")
                  factors <- sapply(model, function(z) myBIC(z,
                    nobs = nobs))
            }
            factors <- factors - min(factors)
            post <- exp(-1/2 * factors)/sum(exp(-1/2 * factors))
        }
        sim.post <- vector("list", length(model))
        f.zeros <- FALSE
        nn <- round(n * post)
        if (verbose)
            print(nn)
        for (i in 1:length(model)) {
            if (nn[i] == 0) {
                f.zeros <- TRUE
                sim.post[[i]] <- sample.qa.posterior(model[[i]],
                  n = 2)
            }
            else {
                sim.post[[i]] <- sample.qa.posterior(model[[i]],
                  n = max(nn[i], 2))
            }
        }
        if (f.zeros == TRUE) {
            warning("One or more models produced zero samples, because of very low posterior probability.")
        }
        if (TRUE) {
            par.names <- sapply(sim.post, function(i) colnames(i))
            upar.names <- unique(unlist(par.names))
            post.avg <- matrix(NA, nrow = sum(nn), ncol = length(upar.names))
            colnames(post.avg) <- upar.names
            current.row <- 1
            for (i in 1:length(model)) {
                if (nn[i] > 0) {
                  start.row <- current.row
                  end.row <- current.row + nn[i] - 1
                  for (j in colnames(sim.post[[i]])) {
                    post.avg[start.row:end.row, j] <- sim.post[[i]][1:nn[i],
                      j]
                  }
                  current.row <- end.row + 1
                }
            }
            post.avg <- as.data.frame(post.avg)
        }
        else {
            post.avg <- sim.post[[1]]
            if (length(model) > 1) {
                for (i in 2:length(model)) {
                  post.avg <- merge(post.avg, sim.post[[i]],
                    all = TRUE, sort = FALSE)
                  if (nn[i] == 0) {
                    nr <- nrow(post.avg)
                    rcut <- (nr - 1):nr
                    post.avg <- post.avg[-rcut, ]
                  }
                }
                if (nn[1] == 0) {
                  post.avg <- post.avg[-(1:2), ]
                }
            }
        }
        if (!is.logical(fill.na))
            post.avg[is.na(post.avg)] <- fill.na
        if (add.names == TRUE) {
            mnames <- match.call()
            mnames <- as.character(mnames[[2]])[2:(length(model) +
                1)]
            rnames <- rep(mnames, times = nn)
            post.avg$model <- rnames
        }
        result <- post.avg
    }
    else {
        mu <- 0
        if (class(model)[1] %in% c("mer", "bmer", "glmerMod",
            "lmerMod")) {
            mu <- fixef(model)
        }
        else {
            mu <- xcoef(model)
        }
        result <- as.data.frame(mvrnorm(n = n, mu = mu, Sigma = vcov(model)))
        if (clean.names == TRUE) {
            for (i in 1:ncol(result)) {
                if (colnames(result)[i] == "(Intercept)") {
                  colnames(result)[i] <- "Intercept"
                }
            }
        }
    }
    result
}

###########################################################################################################
####################################################################### Define Custom Probability Function
dzitnb<-function(x,P1=P1,P2=P2,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {
if (RS[n] == 0){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P2[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P2[n], log = TRUE) + dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}else{

Lp<-increment_log_prob(dbinom(0, size=1, prob=P2[n], log = TRUE) + dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }

###########################################################################################################
####################################################################### Define Custom Probability Function
dthnb<-function(x,P1=P1,P2=P2,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {

Lp<-increment_log_prob(dbinom(ifelse(RS[n]==0,1,0), size=1, prob=P2[n], log = TRUE));

if (RS[n] > 0){
Lp<- increment_log_prob(dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE)  -(1 - log(pnbinom(1, size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = FALSE)))  )  ;
}

}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }


###########################################################################################################
####################################################################### Define Custom Probability Function
dtnb<-function(x,P1=P1,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {
Lp<-increment_log_prob(dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }


################################################################################# Sampling Functions
ZINBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )))*rbinom(1,size=1,prob=(1-logistic(G$Psi2[i] + G$Psi3[i]*scrapExposure[n])));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

#################################################################################### Now Simulate Predictions
ZINBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )))*rbinom(1,size=1,prob=(1-logistic(G$Psi2[i] + G$Psi3[i]*scrapExposure[n])));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

 #################################################################################### Now Simulate Predictions
HNBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
       ################################################################################ Hurdle Model of RS
for (n in 1:N) {
    if(rbinom(1,size=1, prob=logistic(G$Psi2[i] + G$Psi3[i]*(scrapExposure[n])))==1){
          scrapRS[n]<-0} else{

            Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )))
        if(Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapRS[n]<-Nscrap;  }


}

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

#################################################################################### Now Simulate Predictions
HNBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
       ################################################################################ Hurdel Model of RS
for (n in 1:N) {
    if(rbinom(1,size=1, prob=logistic(G$Psi2[i] + G$Psi3[i]*(scrapExposure[n])))==1){
          scrapRS[n]<-0} else{

            Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )))
        if(Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapRS[n]<-Nscrap;  }


}

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }
       
################################################################################# Sampling Functions
NBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }
       
#################################################################################### Now Simulate Predictions
NBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

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

            
####################################################################################### Run Model            
 SkewCalc <- function(RS,Exposure,
MaxExposure=max(Exposure),
Samples=1000,
MethodFit="ETNBShiftingShape",
MethodCI="HPDI",
P=0.95,
Prior="Default",
Updates="TRUE",
PlotResults="TRUE"
){

#############  Priors
if(Prior=="Default"){  
if(MethodFit=="ETZINBShiftingShape" | MethodFit=="ETZINBShiftingShapeScale"){
Priors<<-list( S1=c( 3.0,1.5),   # Scale Param, Log-scale, Constant, Age-Dist Model
               S2=c(-1.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
               S3=c(-1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
               Z1=c( 0.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
               Z2=c(-2.0,1.5),   # Zero inflation Param, Logit-scale, Constant or intecept, RS Model
               Z3=c( 0.0,1.5),   # Zero inflation Param, Logit-scale, slope on exposure, RS Model
               M1=c( 4.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
               M2=c(-2.0,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
               M3=c( 2.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )} else {
if(MethodFit=="ETHNBShiftingShape" | MethodFit=="ETHNBShiftingShapeScale"){
Priors<<-list( S1=c( 3.0,1.5),  # Scale Param, Log-scale, Constant, Age-Dist Model
             S2=c(-1.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
             S3=c( 1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
             Z1=c( 0.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
             Z2=c(-1.0,1.5),   # Zero inflation Param, Logit-scale, Constant or intecept, RS Model
             Z3=c( 0.0,1.5),   # Zero inflation Param, Logit-scale, slope on exposure, RS Model
             M1=c( 5.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
             M2=c( 1.0,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
             M3=c(-1.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )} 
              else {
if(MethodFit=="ETNBShiftingShape" | MethodFit=="ETNBShiftingShapeScale"){
Priors<<-list( S1=c( 4.0,1.5),   # Scale Param, Log-scale, Constant, Age-Dist Model
              S2=c(-2.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
              S3=c(-1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
              Z1=c(-1.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
              M1=c( 5.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
              M2=c( 2.5,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
              M3=c(-1.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )}
             }} }
             else
             {Priors<<-Prior}
########################################################################################################
################################################################################# Define Models
flistShape <- alist(
RS ~ dzitnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

flistShapeScale <- alist(
RS ~ dzitnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

hlistShape <- alist(
RS ~ dthnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

hlistShapeScale <- alist(
RS ~ dthnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

glistShape <- alist(
RS ~ dtnb( P1,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

glistShapeScale <- alist(
RS ~ dtnb( P1,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# END LIBRARY

############################### Begin Function

if(Updates=="TRUE"){
print('Model has been built.')
print('Model is running, please be patient. The parameter estimation process may take several minutes with large datasets.')
}
############################################ Fit Model with MAP
if(MethodFit=="ETZINBShiftingShape"){
fitSkew <- map( flistShape , data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
 } else{
if(MethodFit=="ETZINBShiftingShapeScale"){
 fitSkew <- map(flistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
 }  else {
if(MethodFit=="ETHNBShiftingShape"){
fitSkew <- map( hlistShape , data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETHNBShiftingShapeScale"){
 fitSkew <- map(hlistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETNBShiftingShape"){
 fitSkew <- map(glistShape, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETNBShiftingShapeScale"){
 fitSkew <- map(glistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else(print('Error: Please Select a Valid Method for MethodFit')) }}} }}

if(Updates=="TRUE"){
print('Parameter estimation complete.')
print('Begining sampling from the posterior distribution.')
}

#################################################################################### Now Simulate Predictions
##########################################################################
######################################## Skew Calculation
M_1<-c()
Mc_1<-c()
GGG<-sample.qa.posterior(fitSkew,n=Samples)
FF<-Samples/10
Prog<-round(FF*c(1:10),0)

if(MethodFit=="ETZINBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(ZINBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETZINBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(ZINBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETHNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(HNBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETHNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(HNBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(NBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(NBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
else(print('Error: Please Select a Valid Method for MethodFit')) }}}}}

if(MethodCI=="HPDI"){
  RES<-matrix(NA,ncol=4,nrow=4)

  RES[1,] <-c(M_Mc(c(RS,Exposure))[1],NA,NA,NA)
  RES[2,] <-c(mean(M_1),median( M_1),HPDI(M_1,P)[1],HPDI(M_1,P)[2])

  RES[3,] <-c(M_Mc(c(RS,Exposure))[2],NA,NA,NA)
  RES[4,] <-c(mean(Mc_1),median( Mc_1),HPDI(Mc_1,P)[1],HPDI(Mc_1,P)[2])

  colnames(RES)<-c("MAP","Median","HPDI-L","HPDI-H")
  rownames(RES)<-c("Sample M","Population M","Sample Mc","Population Mc" )

 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning(paste('The Sample Prediction Of M or Mc Is Not Contained in the Models Central' ,(P*100) ,'Percent Posterior Predictive Density.'))
 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning('The Selected Exposure-Truncated Negative-Binomial Model Might Not Be a Good Approximating Model For This Sample.')

} else{
if(MethodCI=="PCI"){

  RES<-matrix(NA,ncol=4,nrow=4)

  RES[1,] <-c(M_Mc(c(RS,Exposure))[1],NA,NA,NA)
  RES[2,] <-c(mean(M_1),median( M_1),PCI(M_1,P)[1],PCI(M_1,P)[2])

  RES[3,] <-c(M_Mc(c(RS,Exposure))[2],NA,NA,NA)
  RES[4,] <-c(mean(Mc_1),median( Mc_1),PCI(Mc_1,P)[1],PCI(Mc_1,P)[2])

  colnames(RES)<-c("MAP","Median","PCI-L","PCI-H")
  rownames(RES)<-c("Sample M","Population M","Sample Mc","Population Mc" )

 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning(paste('The Sample Prediction Of M or Mc Is Not Contained in the Models Central' ,(P*100) ,'Percent Posterior Predictive Density.'))
 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning('The Selected The Exposure-Truncated Zero-Inflated Negative-Binomial Model Might Not Be a Good Approximating Model For This Sample.')

 }
else ( print('Error in Specifying Confidence Interval Method'))
 }
 if(PlotResults=="TRUE"){
print('Plotting requires extra computation, please be patient.')
######################################## Sample Data
library(RColorBrewer)
Samples<-100
FF<-Samples/10
Prog<-round(FF*c(1:10),0)
plotRS<-c()
plotExposure<-c()

if(MethodFit=="ETZINBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-ZINBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
} }

 else{
if(MethodFit=="ETZINBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-ZINBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETHNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-HNBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETHNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-HNBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-NBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-NBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }

else(print('Error: Please Select a Valid Method for MethodFit'))  }}} } }

 
 layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths =c(2,2) , heights=c(1,1))
 plot(density(plotRS,bw=.4), main="RS Densities")
lines(density(RS,bw=.4), col="red")
legend("top", inset=0, title="Distribution",
   c("Sample","Posterior"), fill=c("red","black"), horiz=FALSE)
plot(density(plotExposure,bw=.4), main="Exposure Densities")
lines(density(Exposure,bw=.4), col="red")
legend("top", inset=0, title="Distribution",
   c("Sample","Posterior"), fill=c("red","black"), horiz=FALSE)
   
   

windows()

 ######################################################################## Plot it
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths =c(2,2) , heights=c(1,1))

   jet.colors <-colorRampPalette(c("white",brewer.pal(9,"YlOrRd")))
 MaxRS<-max(c(max(RS,na.rm=T),max(plotRS,na.rm=T)),na.rm=T)
 smoothScatter(RS~Exposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                 colramp = jet.colors,transformation = function(x) x^.45,
                 nbin=200,  xlab="Sample Age Distribution",
                 ylab="Sample RS Distribution")

 smoothScatter(plotRS~plotExposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                colramp = jet.colors,transformation = function(x) x^.45,
                nbin=200, xlab="Population Age Distribution (Predictions)",
                ylab="Predicted Posterior Population RS Distribution")

 }


if(MethodCI=="HPDI"){
print(paste('Confidence Intervals are: ',P*100,'percent Highest Posterior Density Intervals.'))
} else{
if(MethodCI=="PCI"){
print(paste('Confidence Intervals are: Central', P*100,'percent Posterior Credibility Intervals.'))
 }
else ( print('Error in Specifying Confidence Interval Method'))
 }
 print(precis(fitSkew))
 print(RES)
  if(PlotResults=="TRUE"){
return(list(RES=RES,Precis=precis(fitSkew),P_RS=plotRS,P_Exposure=plotExposure))
} else{
return(list(RES=RES,Precis=precis(fitSkew)))
}
 }
 

################################################################################ Example with Fake Data
################ Set Parameters
#N<-100
#MaxExposure<-45

#Kappa1 <- 3
#Kappa2 <- -2
#Psi1 <- -1
#Theta1 <- 5
#Theta2 <- 1.5
#Theta3 <- .5

#scrapExposure<-c()
#scrapRS<-c()

################################################ Split Exposure RNG into point process at MaxExposure and the NB
#Nscrap<-round(logistic(Psi1)*N);
#Npp <- 1;
#while (Npp < Nscrap) {
#Npp <- Npp + 1;
#}
################################################################################ Model of Exposure
#  for( n in 1:Npp){
#       scrapExposure[n]<-MaxExposure;
#       }
#  for( n in (Npp+1):N){
#      Ticker <- 1;
#      while (Ticker == 1) {
#        Nscrap <-rnbinom(1, size=(exp(Theta1)/(exp(Kappa1))), prob=(1 / ( 1 + (exp(Kappa1)) )))
#        if(Nscrap<=MaxExposure){
#        Ticker<-0;}
#        else{
#        Ticker<-1;}}

#       scrapExposure[n]<-Nscrap;
#       }
################################################################################ Model of RS
#   for( n in 1:N){
#   scrapRS[n]<- rnbinom( 1,size=(exp(Theta2+ logistic(Theta3)*log(1+scrapExposure[n]))/ (exp(Kappa2))),prob=(1 / ( 1 + (exp(Kappa2)) )));
#         }
################################################################################ Run SkewCalc
# RES<-SkewCalc(scrapRS,scrapExposure, MethodFit="ETNBShiftingShape")
# Parameters should be recovered
 
################################################################################ Example with Real Data
# data(Kipsigis)
# d<-Kipsigis
# RS<-round(d$soff,0)
# Exposure<-d$age
# Exposure<-round(Exposure-15,0)
# Exposure<-ifelse(Exposure>45,45,Exposure)
# MaxExposure<-max(Exposure)

# RES<-SkewCalc(RS,Exposure,Samples=100, MethodFit="ETNBShiftingShape")
    
 
 
 
