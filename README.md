SkewCalc
========

<img align="right" src="https://github.com/ctross/STRAND/blob/main/logo3.png" alt="logo" width="140"> 

Inequality or skew in reproductive success (RS) is common across many animal species and is of long-standing interest to the study of social evolution. However, the measurement of inequality in RS in natural populations has been challenging because existing quantitative measures are highly sensitive to variation in group/sample size, mean RS, and age-structure. This makes comparisons across multiple groups and/or species vulnerable to statistical artefacts and hinders empirical and theoretical progress. Here, we present a new measure of reproductive skew, the multinomial index, M, that is unaffected by many of the structural biases affecting existing indices. M is analytically related to Nonacs' binomial index, B, and comparably accounts for heterogeneity in age across individuals; in addition, M allows for the possibility of diminishing or even highly nonlinear RS returns to age. Unlike B, however, M is not biased by differences in sample/group size. 
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
?M_index_stan # Explains the printed results

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


4) Finally, an analysis by sex and cultural group:
```{r}
# Load library and attach data
library(ggridges)
library(viridis)
library(SkewCalc) 
data(SukumaMales) 
data(KipsigisMales) 
data(KipsigisFemales) 
data(ColombiaRS) 
d <- ColombiaRS
d$age <- d$age - 13

M_index_stan(d$rs[which(d$group=="AFROCOLOMBIAN" & d$sex=="M")],d$age[which(d$group=="AFROCOLOMBIAN" & d$sex=="M")]) 
M_post_A_male <- extract(StanResults, pars="M_age")$M_age
M_point_A_male <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

M_index_stan(d$rs[which(d$group=="AFROCOLOMBIAN" & d$sex=="F")],d$age[which(d$group=="AFROCOLOMBIAN" & d$sex=="F")]) 
M_post_A_female <- extract(StanResults, pars="M_age")$M_age
M_point_A_female <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

M_index_stan(d$rs[which(d$group=="EMBERA" & d$sex=="M")],d$age[which(d$group=="EMBERA" & d$sex=="M")]) 
M_post_E_male <- extract(StanResults, pars="M_age")$M_age
M_point_E_male <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

M_index_stan(d$rs[which(d$group=="EMBERA" & d$sex=="F")],d$age[which(d$group=="EMBERA" & d$sex=="F")]) 
M_post_E_female <- extract(StanResults, pars="M_age")$M_age
M_point_E_female <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

M_index_stan(KipsigisMales$rs, KipsigisMales$age)
M_post_K_male <- extract(StanResults, pars="M_age")$M_age
M_point_K_male <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

M_index_stan(KipsigisFemales$rs, KipsigisFemales$age)
M_post_K_female <- extract(StanResults, pars="M_age")$M_age
M_point_K_female <- M_index_age(model_dat$r,model_dat$t,Samples=500) 

# Finally, plot the posterior estimates of M by group and sex
df1 <- data.frame(M=M_post_A_male,Sex=rep("Male",length(M_post_A_male)),Group=rep("Afrocolombian",length(M_post_A_male)))
df2 <- data.frame(M=M_post_A_female,Sex=rep("Female",length(M_post_A_female)),Group=rep("Afrocolombian",length(M_post_A_female)))
df3 <- data.frame(M=M_post_E_male,Sex=rep("Male",length(M_post_E_male)),Group=rep("Embera",length(M_post_E_male)))
df4 <- data.frame(M=M_post_E_female,Sex=rep("Female",length(M_post_E_female)),Group=rep("Embera",length(M_post_E_female)))
df6 <- data.frame(M=M_post_K_male,Sex=rep("Male",length(M_post_K_male)),Group=rep("Kipsigis",length(M_post_K_male)))
df7 <- data.frame(M=M_post_K_female,Sex=rep("Female",length(M_post_K_female)),Group=rep("Kipsigis",length(M_post_K_female)))
df <- rbind(df1,df2,df3,df4,df6,df7)

dfb1 <- data.frame(x0=M_point_A_male,Sex=rep("Male",1),Group=rep("Afrocolombian",1))
dfb2 <- data.frame(x0=M_point_A_female,Sex=rep("Female",1),Group=rep("Afrocolombian",1))
dfb3 <- data.frame(x0=M_point_E_male,Sex=rep("Male",1),Group=rep("Embera",1))
dfb4 <- data.frame(x0=M_point_E_female,Sex=rep("Female",1),Group=rep("Embera",1))
dfb6 <- data.frame(x0=M_point_K_male,Sex=rep("Male",1),Group=rep("Kipsigis",1))
dfb7 <- data.frame(x0=M_point_K_female,Sex=rep("Female",1),Group=rep("Kipsigis",1))
dfb <- rbind(dfb1,dfb2,dfb3,dfb4,dfb6,dfb7)

  df$Group <- factor(df$Group)
  df$Group <- factor(df$Group,levels(df$Group)[c(2,1,3)])
  
  dfb$Group <- factor(dfb$Group)
  dfb$Group <- factor(dfb$Group,levels(dfb$Group)[c(2,1,3)])

ggplot()  +
 stat_density_ridges(data=df, aes(x=M, y=Group, fill=0.5 - abs(0.5-..ecdf..)),
  geom = "density_ridges_gradient", calc_ecdf = TRUE, color="white") +
  scale_fill_viridis(name = "Tail probability", direction = -1,option="inferno")+ facet_grid(.~Sex) +   
  theme(strip.text.x = element_text(size=14, face="bold"), strip.text.y = element_text(size=14, face="bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.title=element_text(size=14),legend.text=element_text(size=12))+      
  geom_segment(data=dfb, aes(x = x0, xend = x0, y = as.numeric(Group), yend = as.numeric(Group) + .9), color = "darkred") + 
  theme_ridges(grid = TRUE, center = TRUE) +   geom_hline(yintercept=c(1,2,3),color="white") 


```




