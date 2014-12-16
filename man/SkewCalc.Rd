\name{SkewCalc}
\alias{SkewCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Posterior Confidence Intervals for Estimates of M and Mc with Bayesian Negative Binomial Models}
\description{
SkewCalc allows users to calculate the posterior confidence intervals for estimates of M and Mc with various exposure time truncated Bayesian negative binomial models.
}
\usage{
SkewCalc( RS, Exposure , MaxExposure=max(Exposure), Samples=1000 , MethodFit="ETNBShiftingShape" , MethodCI="HPDI", P=0.95, Prior="Default",Updates="TRUE", PlotResults="TRUE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RS}{A vector of reproductive success (or any count for which skew is to be calculated) values. Must be integer-valued.}
\item{Exposure}{A vector of exposure times to risk of RS. Must be integer-valued.}
\item{MaxExposure}{Upper limit of exposure time. This value allows the model to account for an inflated number of individuals at completed reproduction.}
\item{Samples}{Number of samples from the posterior distribution. Defaults to 1000. }
\item{MethodFit}{Defines the Bayesian model to be used in estimating the population-level joint distribtution of RS and Exposure time. Options are: "ETNBShiftingShape" (Default), "ETNBShiftingShapeScale", "ETZINBShiftingShape", "ETZINBShiftingShapeScale", "ETHNBShiftingShape","ETHNBShiftingShapeScale". See Details.}
\item{MethodCI}{Method for calculating confidences intervals. Options are: "HPDI" (Default) or "PCI". See Details.}
\item{P}{Percentile of confidence interval. Defaults to 0.95.}
\item{Priors}{Priors for the Bayesian model. Defaults depend on the model. See Details.}
\item{Updates}{If "TRUE" prints progress reports.}
\item{PlotResults}{If "TRUE" plots posterior distributions. See Details.}
}
\details{
This command provides a convenient interface for estimating the population-level joint distribtution of RS and Exposure time, using one of a suite of theoretically derived Bayesian models. Posterior estimates of M and Mc are derived from the population-level joint distribtution.
The symbols for \code{MethodFit} stand for "E(xposure) T(runcated) N(egative) B(inomial)", with options for "Z(ero) I(nflated)" and "H(urdle)". "ShiftingShape" implies that the Shape parameter of the Gamma distribution underlying Gamma-Poisson (aka. Negative Binomial) shifts as a function of Exposure time. "ShiftingShapeScale" implies that both the Shape and Scale parameter of the Gamma distribution shifts as a function of Exposure time. 
The symbols for \code{MethodCI} stand for H(ighest) P(osteior) D(ensity) I(nterval) or Central P(osteior) C(redibility) I(nterval).
Custom priors (Gaussian on unconstrained space) can be specified with a list; however, this is not recommended for most users.
If \code{PlotResults}="TRUE", plots of the sample and population-level joint distribution, and density plots of the sample and population-level RS and Exposure variables, will be returned with the estimates of M and Mc. These plots should appear similar. Divergence between the sample and poulation estimates might be indicative that a better model is needed.
}
\value{
Returns a list with following items:
\item{RES}{Sample and population estimates of M and Mc.}
\item{Precis}{The estimated posterior densities of the model parameters.}
\item{P_RS}{If \code{PlotResults}="TRUE", estimates of RS are returned so that the user can use their own plotting methods if desired.}
\item{P_Exposure}{If \code{PlotResults}="TRUE", estimates of Exposure are returned so that the user can use their own plotting methods if desired.}

}
\references{Hooper, Ross, etc... ProcB?}
\author{Cody T. Ross}
\seealso{\code{Rethinking} by Richard McElreath for Maximum A Posteriori methods.}
\examples{
################################################################################ Example with Fake Data
################ Set Parameters
N<-100
MaxExposure<-45

Kappa1 <- 5
Kappa2 <- .5
Psi1 <- -1
Theta1 <- 5
Theta2 <- 1.5
Theta3 <- .5

scrapExposure<-c()
scrapRS<-c()

################################################ Split Exposure RNG into point process at MaxExposure and the NB
Nscrap<-round(logistic(Psi1)*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}
################################################################################ Model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(Theta1)/(exp(Kappa1))), prob=(1 / ( 1 + (exp(Kappa1)) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }
################################################################################ Model of RS
   for( n in 1:N){
scrapRS[n]<- rnbinom( 1,size=(exp(Theta2+ logistic(Theta3)*log(1+scrapExposure[n]))/ (exp(Kappa2))),prob=(1 / ( 1 + (exp(Kappa2)) )));
         }
################################################################################ Run SkewCalc
 RES<-SkewCalc(scrapRS,scrapExposure, MethodFit="ETNBShiftingShape")

# Parameters should be recovered

 
################################################################################ Example with Real Data
 \dontrun{data(Kipsigis)
 d<-Kipsigis
 RS<-round(d$soff,0)
 Exposure<-d$age
 Exposure<-round(Exposure-15,0)
 Exposure<-ifelse(Exposure>45,45,Exposure)
 MaxExposure<-max(Exposure)

 RES<-SkewCalc(RS,Exposure,Samples=100, MethodFit="ETNBShiftingShape")}
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
