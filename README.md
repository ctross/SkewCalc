SkewCalc
========
This is an R package for calculating reproductive skew.

Install by running on R:
```{r}
library(devtools)
install_github("Ctross/SkewCalc")
```

Some quick examples:

1) We use a small data-set from Monique Borgerhoff Mulder collected from Sukuma Men. 
```{r}
# Load library and attach data
library(SkewCalc)  
data(SukumaMales) 

# Prepare data
RS <- SukumaMales$rs  # RS Data 
Time <- SukumaMales$age - 18  # Exposure time data. Age mininus minimum age of reproduction. 

# Fit models
M_index_stan(RS,Time) 
?M_index_stan # To explains the printed results

# Contrast the Stan model to the point estimates
M_index(RS,Time) 
M_index_age(RS,Time) 
M_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

Mraw_index(RS,Time) 
Mraw_index_age(RS,Time) 
Mraw_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

# Finally, plot the posterior estimates of M or Mraw
skew_index_plot("M",Age=FALSE)
skew_index_plot("M",Age=TRUE)
skew_index_plot("Mraw",Age=FALSE)
skew_index_plot("Mraw",Age=TRUE)
```

2) We use a larger data-set from Monique Borgerhoff Mulder collected from Kipsigis Men. 
```{r}
# Load library and attach data
library(SkewCalc)  
data(KipsigisMales) 

# Prepare data
RS <- KipsigisMales$rs  # RS Data 
Time <- KipsigisMales$age - min(KipsigisMales$age) + 1  # Exposure time data. Age mininus minimum age of reproduction, plus 1. 

# Fit models
M_index_stan(RS,Time) 

# Contrast the Stan model to the point estimates
M_index(RS,Time) 
M_index_age(RS,Time) 
M_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

Mraw_index(RS,Time) 
Mraw_index_age(RS,Time) 
Mraw_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

# Finally, plot the posterior estimates of M or Mraw
skew_index_plot("M",Age=FALSE)
skew_index_plot("M",Age=TRUE)
skew_index_plot("Mraw",Age=FALSE)
skew_index_plot("Mraw",Age=TRUE)
```


3) Now, lets make up some fake data and test the models.
```{r}
# Load library and attach data
library(SkewCalc)  

# Prepare data
N <- 1000
Alpha <- 1.1
Gamma <- 0.81

Time <- round(runif(N,1,70),0)

RS <- rep(NA,N)
for(i in 1:N)
RS[i] <- rpois(1,Alpha*Time[i]^Gamma)

# Fit models
M_index_stan(RS,Time) 

# Contrast the Stan model to the point estimates
M_index(RS,Time) 
M_index_age(RS,Time) 
M_index_from_B_index(B_index(RS,Time),sum(RS),length(RS))

Mraw_index(RS,Time) 
Mraw_index_age(RS,Time) 
Mraw_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

# Finally, plot the posterior estimates of M or Mraw
skew_index_plot("M",Age=FALSE)
skew_index_plot("M",Age=TRUE)
skew_index_plot("Mraw",Age=FALSE)
skew_index_plot("Mraw",Age=TRUE)
```

4) Now, lets introduce more skew.
```{r}
# Load library and attach data
library(SkewCalc)  

# Prepare data
N <- 1000
Alpha <- 1.1
Gamma <- 0.81

Time <- round(runif(N,1,70),0)

RS <- rep(NA,N)
for(i in 1:N){
Mu <- Alpha*Time[i]^Gamma
B <- 0.1
Lambda <- rgamma(1,Mu*B,B)
RS[i] <- rpois(1,Lambda)
}

# Fit models
M_index_stan(RS,Time) 

# Contrast the Stan model to the point estimates
M_index(RS,Time) 
M_index_age(RS,Time) 
M_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

Mraw_index(RS,Time) 
Mraw_index_age(RS,Time) 
Mraw_index_from_B_index(B_index(RS,Time),sum(RS),length(RS)) 

# Finally, plot the posterior estimates of M or Mraw
skew_index_plot("M",Age=FALSE)
skew_index_plot("M",Age=TRUE)
skew_index_plot("Mraw",Age=FALSE)
skew_index_plot("Mraw",Age=TRUE)
```







