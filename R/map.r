#This code is lifted directly from Richard McElreath, Rethninking https://github.com/rmcelreath/rethinking
# Create Function Here, Create Class Below
map<-function (flist, data, start, method = "Nelder-Mead", hessian = TRUE,
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
