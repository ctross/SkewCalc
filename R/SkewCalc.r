#
# If you don't have stan, download and install R-Tools 
# https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows
#
# Then follow instructions https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started, or just run
# source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
# install_rstan()
#

################################################################################
##### Begin Model Code in Stan
################################################################################
skew_code_GRF<-'
################################################################################
##### Declare functions to be used in model
################################################################################
functions{
############################### Function to calculate M (slot 1) and Mc (slot 2)
 vector M(vector Scrap, int N) {
 real K;                           #  Declare types
 real F;                           #
 vector[N] Si;                     #
 real S;                           #
 vector[2] M;                      #
 vector[N] scrapRS;                #
 vector[N] scrapExposure;          #

  for(n in 1:N){                   # Split input into RS and Exposure vectors
   scrapRS[n]<-Scrap[n];           #
   scrapExposure[n]<-Scrap[N+n];   #
   }                               #

############################### Calculate M and Mc
	K <- sum(scrapRS);
	F <- sum(scrapExposure);
	Si <- ((scrapRS/K)-(scrapExposure/F)) .* ((scrapRS/K)-(scrapExposure/F));
	S <- sum(Si);
 	M[1] <- sqrt(N * S);                      # M
	M[2] <- sqrt(N * S)*sqrt(mean(scrapRS));  # Mc
  return M;
}

############################### Function to simulate predictions from the model
vector M_NB_rng(int N, int MaxExposure, vector Theta, vector GRF_RS_Mu, vector GRF_RS_Zero){
 int scrapRS[N];            # Declare types
 int scrapExposure[N];      #
 vector[N*2] scrap;         #
                            #
 real NscrapR;              #
 int Npp;                   #
 int Ticker;                #
                            #
 int Nscrap;                #
 real P1;                   #
 vector[N] Mu1;             #
 vector[N] B1;              #


############################
P1 <- inv_logit(Theta[1]);        # Model of Point Density on Exposure=max(Exposure)
for(n in 1:N){
Mu1[n]<- exp(Theta[2]);           # Model of Mean of Exposure Model
B1[n] <- exp(Theta[3]);           # Model of Scale of Exposure Model
      }

############################ Split RNG into point process at MaxExp and the ZINB
NscrapR<-round(P1*N);
Npp <- 1;
while (Npp <= NscrapR) {
Npp <- Npp + 1;
}
# Npp is the number of slots allocated to the density on max(Exposure) by P1

############################################################################ RNG
################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
       Ticker <- 1;
       while (Ticker == 1) {
        Nscrap <-neg_binomial_rng(Mu1[n]*B1[n],B1[n]);
        if(Nscrap<=MaxExposure && Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}
        }
        scrapExposure[n]<-Nscrap;
       }
##################################################### RS conditional on exposure
  for( n in 1:N){
    scrapRS[n]<-neg_binomial_rng(exp(Theta[4] + GRF_RS_Mu[scrapExposure[n]])*exp(Theta[6]),exp(Theta[6]));
       }
#################################################### Store samples in one vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n]*bernoulli_rng((1-inv_logit(Theta[7] + GRF_RS_Zero[scrapExposure[n]])));
       scrap[N+n]<-scrapExposure[n];
       }
  return scrap;
       }

}

################################################################################
##### Declare data to be used in model
################################################################################
data {
int<lower=0> N;                # Number sampled individuals
int<lower=0> RS[N];            # Reproductive success
int<lower=0> Exposure[N];      # Time at risk for RS
int<lower=0> MaxExposure;      # Max time at risk for RS
vector[MaxExposure] ZeroPrior; # GFR prior mean
}

################################################################################
##### Transform some data to vectors to speed up downstream calculations
################################################################################
transformed data {
vector<lower=0>[N] RS_v;            # Reproductive Success
vector<lower=0>[N] Exposure_v;      # Time at Risk for RS

for(n in 1:N){
 RS_v[n]<-RS[n];
 Exposure_v[n]<-Exposure[n];
             }           
}

################################################################################
##### Declare parameters to be estimated in the model
################################################################################
parameters {
vector[12] Theta;
vector[MaxExposure] GRF_RS_Mu;
vector[MaxExposure] GRF_RS_Zero;
}

################################################################################
##### Fit model
################################################################################
model {
################### Declare local storage arrays
real G;
real ME;
matrix[MaxExposure,MaxExposure] X;

vector[N] Mu1;
vector[N] Mu2;

vector[N] B1;
vector[N] B2;

real P1;
vector[N] P2;

######################################################################### Priors
Theta ~ normal(0,5);

######################################################### Gaussian Random Fields
####### Construct Gaussian Process for Mu
for (i in 1:(MaxExposure-1)){
for (j in (i+1):MaxExposure){
 ME <- MaxExposure;
  G <-((j-i) * (j-i))/(ME*ME);                                                  # Calculate distance
                X[i,j] <- inv_logit(Theta[9]) * exp( -exp(Theta[10]) * G);      # Estimate Correlations
                X[j,i] <- X[i,j];                                               # Fill Other Triangle
                       }}

for (i in 1:MaxExposure){
                X[i,i] <- 1;                                               # Fill Diag
                   }

         X <- cholesky_decompose(X);                                            # Decompose Rho

         GRF_RS_Mu ~ multi_normal_cholesky( ZeroPrior , exp(Theta[5])*X );      # Model GP in Cholesky notation for efficiancy

####### Construct Gaussian Process for Zeros
for (i in 1:(MaxExposure-1)){
for (j in (i+1):MaxExposure){
 ME <- MaxExposure;
  G <-((j-i) * (j-i))/(ME*ME);                                                  # Calculate distance
                X[i,j] <- inv_logit(Theta[11]) * exp( -exp(Theta[12]) * G);     # Estimate Correlations
                X[j,i] <- X[i,j];                                               # Fill Other Triangle
                       }}

for (i in 1:MaxExposure){
                X[i,i] <- 1;                                                    # Fill Diag
                   }

         X <- cholesky_decompose(X);                                            # Decompose Rho

         GRF_RS_Zero ~ multi_normal_cholesky( ZeroPrior , exp(Theta[8])*X );    # Model GP in Cholesky notation for efficiancy

############################## Model of Point Density on Exposure=max(Exposure)
P1 <- inv_logit(Theta[1]);

################################################################################
for(n in 1:N){ 
############################################## Link Functions for Exposure Model
Mu1[n]<- exp(Theta[2]);           # Model of Mean of Exposure Model
B1[n] <- exp(Theta[3]);           # Model of Scale of Exposure Model

##################################################### Link Function for RS Model
Mu2[n]<- exp(Theta[4] + GRF_RS_Mu[Exposure[n]]);                  # Model of Mean of RS Model
B2[n] <- exp(Theta[6] );                                          # Model of Scale of RS Model
P2[n] <- inv_logit(Theta[7] + GRF_RS_Zero[Exposure[n]]);          # Model of Zeros of RS Model
}

############################################ Negative Binomial Model of Exposure
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
increment_log_prob( log_sum_exp(bernoulli_log(1,P1), bernoulli_log(0,P1) + (neg_binomial_log(Exposure[n], Mu1[n]*B1[n], B1[n])-neg_binomial_cdf_log(MaxExposure, Mu1[n]*B1[n], B1[n])))); 
}else {
increment_log_prob( bernoulli_log(0,P1) + (neg_binomial_log(Exposure[n], Mu1[n]*B1[n],  B1[n])-neg_binomial_cdf_log(MaxExposure, Mu1[n]*B1[n], B1[n]))); 
}
}
################################################## Negative Binomial Model of RS
for (n in 1:N) {
if(RS[n]==0){
increment_log_prob( log_sum_exp(bernoulli_log(1,P2[n]), bernoulli_log(0,P2[n]) + neg_binomial_log(RS[n], Mu2[n]*B2[n],  B2[n]) )); 
}else {
increment_log_prob( bernoulli_log(0,P2[n]) + neg_binomial_log(RS[n], Mu2[n]*B2[n],  B2[n])); 
}
}
}

################################################################################
##### Simulate from the posterior
################################################################################
generated quantities{
vector[2] M_Mc;
vector[2*N] Pred;

Pred <- M_NB_rng(N,  MaxExposure, Theta, GRF_RS_Mu, GRF_RS_Zero);
M_Mc  <- M(Pred ,N);
}
'

skew_code_Fast<-'
################################################################################
##### Declare functions to be used in model
################################################################################
functions{
############################### Function to calculate M (slot 1) and Mc (slot 2)
 vector M(vector Scrap, int N) {
 real K;                           #  Declare types
 real F;                           #
 vector[N] Si;                     #
 real S;                           #
 vector[2] M;                      #
 vector[N] scrapRS;                #
 vector[N] scrapExposure;          #

  for(n in 1:N){                   # Split input into RS and Exposure vectors
   scrapRS[n]<-Scrap[n];           #
   scrapExposure[n]<-Scrap[N+n];   #
   }                               #

############################### Calculate M and Mc
	K <- sum(scrapRS);
	F <- sum(scrapExposure);
	Si <- ((scrapRS/K)-(scrapExposure/F)) .* ((scrapRS/K)-(scrapExposure/F));
	S <- sum(Si);
 	M[1] <- sqrt(N * S);                      # M
	M[2] <- sqrt(N * S)*sqrt(mean(scrapRS));  # Mc
  return M;
}

############################### Function to simulate predictions from the model
vector M_NB_rng(int N, int MaxExposure, vector Theta){
 int scrapRS[N];            # Declare types
 int scrapExposure[N];      #
 vector[N*2] scrap;         #
                            #
 real NscrapR;              #
 int Npp;                   #
 int Ticker;                #
                            #
 int Nscrap;                #
 real P1;                   #
 vector[N] Mu1;             #
 vector[N] B1;              #


############################
P1 <- inv_logit(Theta[1]);        # Model of Point Density on Exposure=max(Exposure)
for(n in 1:N){
Mu1[n]<- exp(Theta[2]);           # Model of Mean of Exposure Model
B1[n] <- exp(Theta[3]);           # Model of Scale of Exposure Model
      }

############################ Split RNG into point process at MaxExp and the ZINB
NscrapR<-round(P1*N);
Npp <- 1;
while (Npp <= NscrapR) {
Npp <- Npp + 1;
}
# Npp is the number of slots allocated to the density on max(Exposure) by P1

############################################################################ RNG
################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
       Ticker <- 1;
       while (Ticker == 1) {
        Nscrap <-neg_binomial_rng(Mu1[n]*B1[n],B1[n]);
        if(Nscrap<=MaxExposure && Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}
        }
        scrapExposure[n]<-Nscrap;
       }
##################################################### RS conditional on exposure
  for( n in 1:N){
    scrapRS[n]<-neg_binomial_rng(exp(Theta[4] + Theta[5]*log(scrapExposure[n]))*exp(Theta[6]),exp(Theta[6]));
       }
#################################################### Store samples in one vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n]*bernoulli_rng((1-inv_logit(Theta[7] + Theta[8]*scrapExposure[n])));
       scrap[N+n]<-scrapExposure[n];
       }
  return scrap;
       }

}

################################################################################
##### Declare data to be used in model
################################################################################
data {
int<lower=0> N;                # Number sampled individuals
int<lower=0> RS[N];            # Reproductive success
int<lower=0> Exposure[N];      # Time at risk for RS
int<lower=0> MaxExposure;      # Max time at risk for RS
vector[MaxExposure] ZeroPrior; # GFR prior mean
}

################################################################################
##### Transform some data to vectors to speed up downstream calculations
################################################################################
transformed data {
vector<lower=0>[N] RS_v;            # Reproductive Success
vector<lower=0>[N] Exposure_v;      # Time at Risk for RS

for(n in 1:N){
 RS_v[n]<-RS[n];
 Exposure_v[n]<-Exposure[n];
             }
}

################################################################################
##### Declare parameters to be estimated in the model
################################################################################
parameters {
vector[8] Theta;
}

################################################################################
##### Fit model
################################################################################
model {
################### Declare local storage arrays
real G;
real ME;
matrix[MaxExposure,MaxExposure] X;

vector[N] Mu1;
vector[N] Mu2;

vector[N] B1;
vector[N] B2;

real P1;
vector[N] P2;

######################################################################### Priors
Theta ~ normal(0,5);

############################## Model of Point Density on Exposure=max(Exposure)
P1 <- inv_logit(Theta[1]);

################################################################################
for(n in 1:N){
############################################## Link Functions for Exposure Model
Mu1[n]<- exp(Theta[2]);           # Model of Mean of Exposure Model
B1[n] <- exp(Theta[3]);           # Model of Scale of Exposure Model

##################################################### Link Function for RS Model
Mu2[n]<- exp(Theta[4] + Theta[5]*log(Exposure[n]));               # Model of Mean of RS Model
B2[n] <- exp(Theta[6] );                                          # Model of Scale of RS Model
P2[n] <- inv_logit(Theta[7] + Theta[8]*Exposure[n]);              # Model of Zeros of RS Model
}

############################################ Negative Binomial Model of Exposure
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
increment_log_prob( log_sum_exp(bernoulli_log(1,P1), bernoulli_log(0,P1) + (neg_binomial_log(Exposure[n], Mu1[n]*B1[n], B1[n])-neg_binomial_cdf_log(MaxExposure, Mu1[n]*B1[n], B1[n]))));
}else {
increment_log_prob( bernoulli_log(0,P1) + (neg_binomial_log(Exposure[n], Mu1[n]*B1[n],  B1[n])-neg_binomial_cdf_log(MaxExposure, Mu1[n]*B1[n], B1[n])));
}
}
################################################## Negative Binomial Model of RS
for (n in 1:N) {
if(RS[n]==0){
increment_log_prob( log_sum_exp(bernoulli_log(1,P2[n]), bernoulli_log(0,P2[n]) + neg_binomial_log(RS[n], Mu2[n]*B2[n],  B2[n]) ));
}else {
increment_log_prob( bernoulli_log(0,P2[n]) + neg_binomial_log(RS[n], Mu2[n]*B2[n],  B2[n]));
}
}
}

################################################################################
##### Simulate from the posterior
################################################################################
generated quantities{
vector[2] M_Mc;
vector[2*N] Pred;

Pred <- M_NB_rng(N,  MaxExposure, Theta);
M_Mc  <- M(Pred ,N);
}
'

M_index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	C <- sqrt(N * S)
  C
}

Mc_index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	D1 <- sqrt(N * S)*sqrt(mean(ki))
  D1
}

 concat<-function (...)
{
    paste(..., collapse = "", sep = "")
}

HPDI<-function (samples, prob = 0.89)
{
    class.samples <- class(samples)[1]
    coerce.list <- c("numeric", "matrix", "data.frame", "integer",
        "array")
    if (class.samples %in% coerce.list) {
        samples <- coda::as.mcmc(samples)
    }
    x <- sapply(prob, function(p) coda::HPDinterval(samples,
        prob = p))
    n <- length(prob)
    result <- rep(0, n * 2)
    for (i in 1:n) {
        low_idx <- n + 1 - i
        up_idx <- n + i
        result[low_idx] <- x[1, i]
        result[up_idx] <- x[2, i]
        names(result)[low_idx] <- concat("|", prob[i])
        names(result)[up_idx] <- concat(prob[i], "|")
    }
    return(result)
}


SkewCalc<-function(RS,Exposure, Samples=1000, Warmup=500, Chains=1, Refresh=1, Code="GRF"){
################## Load libraries
 require(rstan)

################## Clean inputs
 SCRAP<-data.frame(RS,Exposure)
 SCRAP<-SCRAP[complete.cases(SCRAP),]

 RS<-SCRAP$RS
 Exposure<-SCRAP$Exposure

 N <- length(RS)
 MaxExposure<-max(Exposure)
 MaxRS<-max(RS)
 ZeroPrior<-rep(0, MaxExposure)

################## Prepare data for Stan
 model_dat<-list(
  N=N,
  RS=RS,
  Exposure=Exposure,
  MaxExposure=MaxExposure,
  ZeroPrior=ZeroPrior
  )

  if(Code=="GRF"){
 StanResults <- stan(model_code=skew_code_GRF, data=model_dat, thin=1, iter=Samples, warmup=Warmup, chains=Chains, refresh=Refresh)
 MMC<-extract(StanResults, pars="M_Mc")$M_Mc

 SkewFit<-matrix(NA,nrow=4,ncol=3)
 SkewFit[1,1]<-M_index(RS,Exposure)
 SkewFit[2,1]<-Mc_index(RS,Exposure)
 SkewFit[3,1]<-median(MMC[,1])
 SkewFit[4,1]<-median(MMC[,2])
 SkewFit[3,2]<-HPDI(MMC[,1])[1]
 SkewFit[4,2]<-HPDI(MMC[,2])[1]
 SkewFit[3,3]<-HPDI(MMC[,1])[2]
 SkewFit[4,3]<-HPDI(MMC[,2])[2]

 colnames(SkewFit)<-c("Estimate","Lower_90%_HPDI","Upper_90%_HPDI")
 rownames(SkewFit)<-c("M_Point_Estimate","Mc_Point_Estimate","M_Posterior_Estimate","Mc_Posterior_Estimate")

 print(SkewFit)
 SkewResults<<-list(SkewFit=SkewFit,StanResults=StanResults)}
 else{   if(Code=="Fast"){
 StanResults <- stan(model_code=skew_code_Fast, data=model_dat, thin=1, iter=Samples, warmup=Warmup, chains=Chains, refresh=Refresh)
       MMC<-extract(StanResults, pars="M_Mc")$M_Mc

 SkewFit<-matrix(NA,nrow=4,ncol=3)
 SkewFit[1,1]<-M_index(RS,Exposure)
 SkewFit[2,1]<-Mc_index(RS,Exposure)
 SkewFit[3,1]<-median(MMC[,1])
 SkewFit[4,1]<-median(MMC[,2])
 SkewFit[3,2]<-HPDI(MMC[,1])[1]
 SkewFit[4,2]<-HPDI(MMC[,2])[1]
 SkewFit[3,3]<-HPDI(MMC[,1])[2]
 SkewFit[4,3]<-HPDI(MMC[,2])[2]

 colnames(SkewFit)<-c("Estimate","Lower_90%_HPDI","Upper_90%_HPDI")
 rownames(SkewFit)<-c("M_Point_Estimate","Mc_Point_Estimate","M_Posterior_Estimate","Mc_Posterior_Estimate")

 print(SkewFit)
 SkewResults<<-list(SkewFit=SkewFit,StanResults=StanResults)}
       else{
       testit <- function() warning("Code variable not defined correctly", call. = FALSE)
        testit()
}}
        }


skew_diagnositic_plots<-function(RS,Exposure,SkewResults=SkewResults){
 require(RColorBrewer)
 Pred<-extract(SkewResults$StanResults, pars="Pred")$Pred
 N<-length(RS)
 MaxExposure<-max(Exposure)
 plotRS <- c(Pred[,1:N])
 plotExposure <- c(Pred[,(1+N):(2*N)])

 layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths =c(2,2) , heights=c(1,1))
 plot(density(plotRS,bw=2), main="RS Densities")
 lines(density(RS,bw=2), col="red")
 legend("top", inset=0, title="Distribution",
 c("Sample","Posterior"), fill=c("red","black"), horiz=FALSE)
 plot(density(plotExposure,bw=2), main="Exposure Densities")
 lines(density(Exposure,bw=2), col="red")
 legend("top", inset=0, title="Distribution",
 c("Sample","Posterior"), fill=c("red","black"), horiz=FALSE)

 windows()

 layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths =c(2,2) , heights=c(1,1))
 jet.colors <-colorRampPalette(c("white",brewer.pal(9,"YlOrRd")))
 MaxRS<-max(c(max(RS,na.rm=T),max(plotRS,na.rm=T)),na.rm=T)
 smoothScatter(RS~Exposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                 colramp = jet.colors,transformation = function(x) x^.3,  bandwidth=c(1,1),
                 nbin=200,  xlab="Sample Age Distribution",
                 ylab="Sample RS Distribution")
                 points(RS~Exposure,pch=18)
 smoothScatter(plotRS~plotExposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                colramp = jet.colors,transformation = function(x) x^.3, bandwidth=c(1,1),
                nbin=200, nrpoints =0,xlab="Population Age Distribution (Predictions)",
                ylab="Predicted Posterior Population RS Distribution")
                 points(RS~Exposure,pch=18)
}


skew_posterior_plot<-function(RS,Exposure,SkewResults=SkewResults){
 require(RColorBrewer)
  jet.colors <-colorRampPalette(c("white",brewer.pal(9,"YlOrRd")))
 Pred<-extract(SkewResults$StanResults, pars="Pred")$Pred
  N<-length(RS)
  MaxExposure<-max(Exposure)
 plotRS <- c(Pred[,1:N])
 plotExposure <- c(Pred[,(1+N):(2*N)])
  MaxRS<-max(c(max(RS,na.rm=T),max(plotRS,na.rm=T)),na.rm=T)

 smoothScatter(plotRS~plotExposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                colramp = jet.colors,transformation = function(x) x^.3, bandwidth=c(1,1),
                nbin=200, nrpoints =0,xlab="Population Age Distribution (Predicted)",
                ylab="Posterior Population RS Distribution (Predicted)")
                points(RS~Exposure,pch=18)
}

