\name{skew_diagnositic_plots}
\alias{skew_diagnositic_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnositic Plots of Model Fit from Skew Calc}
\description{
To be run after SkewCalc. Displays two sets of plot. The first set plots the univariate RS and Exposure time densities, both the sample densities and the posterior densities. The second set plots the sample and posterior bivariate relationship between RS and Exposure time.
}
\usage{
skew_diagnositic_plots(RS,Exposure,SkewResults=SkewResults)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RS}{A vector of reproductive success (or any count for which skew is to be calculated) values. Must be integer-valued.}
\item{Exposure}{A vector of exposure times to risk of RS. Must be integer-valued.}
\item{SkewResults}{Generated automatically from SkewCalc. }
}
\details{
}
\value{
Two sets of plots are displayed.
}
\references{}
\author{Cody Ross}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
