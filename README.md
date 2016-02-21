SkewCalc
========
This is an R package for Calculating reproductive skew.

Install by running on R:
```{r}
library(devtools)
install_github("Ctross/SkewCalc")
```

A quick example:
```{r}
library(SkewCalc)
SkewCalc(SampleRS,SampleExposure, Samples=400, Warmup=200, Chains=1, Refresh=1, Code="Fast")
```
SkewResults$SkewFit
SkewResults$StanResults
skew_diagnositic_plots(SampleRS,SampleExposure,SkewResults)




