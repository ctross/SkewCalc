
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

############################### Calculate Mraw and M
	K <- sum(scrapRS);
	F <- sum(scrapExposure);
	if(K==0 || F==0){
	M[1]<-99999999;
	M[2]<-99999999;
 return M;	
	}else{
	Si <- ((scrapRS/K)-(scrapExposure/F)) .* ((scrapRS/K)-(scrapExposure/F));
	S <- sum(Si);
 	M[1] <- sqrt(N * S);                      # M
	M[2] <- sqrt(N * S)*sqrt(mean(scrapRS));  # Mc
  return M;
  }
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

if(Npp >= N){
Npp<-N;

  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
       }else{
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
Theta ~ normal(0,3);

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
vector[2] Mraw_M;
vector[2*N] Pred;

Pred <- M_NB_rng(N,  MaxExposure, Theta, GRF_RS_Mu, GRF_RS_Zero);
Mraw_M  <- M(Pred ,N);
}
'
