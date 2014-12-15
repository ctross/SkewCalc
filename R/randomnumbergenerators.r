# These RNGs are written by Cody T Ross and correspond to the density functions used in the package
################################################################################# Sampling Functions
ZINBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )))*rbinom(1,size=1,prob=(1-logistic(G$Psi2[i] + G$Psi3[i]*scrapExposure[n])));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

#################################################################################### Now Simulate Predictions
ZINBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )))*rbinom(1,size=1,prob=(1-logistic(G$Psi2[i] + G$Psi3[i]*scrapExposure[n])));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

 #################################################################################### Now Simulate Predictions
HNBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
       ################################################################################ Hurdle Model of RS
for (n in 1:N) {
    if(rbinom(1,size=1, prob=logistic(G$Psi2[i] + G$Psi3[i]*(scrapExposure[n])))==1){
          scrapRS[n]<-0} else{

            Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )))
        if(Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapRS[n]<-Nscrap;  }


}

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }

#################################################################################### Now Simulate Predictions
HNBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
       ################################################################################ Hurdel Model of RS
for (n in 1:N) {
    if(rbinom(1,size=1, prob=logistic(G$Psi2[i] + G$Psi3[i]*(scrapExposure[n])))==1){
          scrapRS[n]<-0} else{

            Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )))
        if(Nscrap>0){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapRS[n]<-Nscrap;  }


}

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }
       
################################################################################# Sampling Functions
NBShape_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp(G$Kappa2[i]))),prob=(1 / ( 1 + (exp(G$Kappa2[i])) )));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }
       
#################################################################################### Now Simulate Predictions
NBShapeScale_rng<-function(i, G=sample.qa.posterior(fitSkew,n=Samples), MaxExposure, N=length(RS)) {
##################################### Declarations

scrapExposure<-c()
scrapRS<-c()
scrap<-c()
################################################ Split RNG into point process at MaxExp and the ZINB
Nscrap<-round(logistic(G$Psi1[i])*N);
Npp <- 1;
while (Npp < Nscrap) {
Npp <- Npp + 1;
}

########################################################################### RNG
################################################################################ Inflation model of Exposure
  for( n in 1:Npp){
       scrapExposure[n]<-MaxExposure;
       }
  for( n in (Npp+1):N){
      Ticker <- 1;
      while (Ticker == 1) {
        Nscrap <-rnbinom(1, size=(exp(G$Theta1[i])/(exp(G$Kappa1[i]))), prob=(1 / ( 1 + (exp(G$Kappa1[i])) )))
        if(Nscrap<=MaxExposure){
        Ticker<-0;}
        else{
        Ticker<-1;}}

       scrapExposure[n]<-Nscrap;
       }

################################################################################ RS Conditional on Exposure
  for( n in 1:N){
       scrapRS[n]<-rnbinom( 1,size=(exp(G$Theta2[i]+ logistic(G$Theta3[i])*log(1+scrapExposure[n]))/ (exp( G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]) ))),prob=(1 / ( 1 + (exp(G$Kappa2[i] + G$Kappa3[i]*log(1+scrapExposure[n]))) )));

       }

################################################################################ Store Samples in One Vector
  for(n in 1:N){
       scrap[n]<-scrapRS[n];
       scrap[N+n]<-scrapExposure[n];
       }

  scrap

       }
