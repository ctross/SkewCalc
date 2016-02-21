SkewCalc
========
This is an R package for calculating reproductive skew.

Install by running on R:
```{r}
library(devtools)
install_github("Ctross/SkewCalc")
```

A quick example:
```{r}
library(SkewCalc)
SkewCalc(SampleRS,SampleExposure, Samples=400, Warmup=200, Chains=1, Refresh=1, Code="Fast")

# Print results
SkewResults$SkewFit

# Print MCMC summary stats and convergence diagnostics
SkewResults$StanResults

# Plot diagnostics of model fit
skew_diagnositic_plots(SampleRS,SampleExposure,SkewResults)
```




