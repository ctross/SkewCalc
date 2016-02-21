\name{SkewCalc}
\alias{SkewCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Posterior Confidence Intervals for Estimates of M and Mc with Bayesian Negative Binomial Models}
\description{
SkewCalc allows users to calculate the posterior confidence intervals for estimates of M and Mc with a Bayesian exposure-time truncated negative binomial model.
}
\usage{
SkewCalc(RS,Exposure, Samples=1000, Warmup=500, Chains=1, Refresh=1, Code="GRF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RS}{A vector of reproductive success (or any count for which skew is to be calculated) values. Must be integer-valued.}
\item{Exposure}{A vector of exposure times to risk of RS. Must be integer-valued.}
\item{Samples}{Number of samples from the posterior distribution. Defaults to 1000. }
\item{Warmup}{Number of Samples used for warmup. Defaults to 500. Must be < Samples. }
\item{Chain}{Number of MCMC chains to run. Defaults to 1. }
\item{Refresh}{Refresh rate for MCMC progress, defaults to 1.}
\item{Code}{Defaults to "GRF", or a Gassian Random Field model. Can be set to "Fast", for a quicker running model.}
}
\details{
This command provides a convenient interface for estimating the population-level joint distribtution of RS and Exposure time, using a theoretically derived Bayesian model. See https://github.com/Ctross/SkewCalc for model details.
install.packages("http://cran.r-project.org/src/contrib/Archive/RNetLogo/RNetLogo_0.9-6.tar.gz", repo=NULL, type="source")}
\value{
Returns a list with following items:
\item{SkewFit}{Sample and population estimates of M and Mc.}
\item{StanResults}{The Stanfit object from MCMC sampling in Rstan.}
}
\references{Hooper, Ross, etc... ProcB?}
\author{Cody T. Ross}
\examples{
################################################################################ Example with Fake Data
}
\keyword{ }
