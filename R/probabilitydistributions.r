# These functions are custom probability densities, written by Cody T. Ross

####################################################################### Define Custom Probability Function
#Exposure Truncated Negative Binomial, Zero Inflated RS Negative Binomial
dzitnb<-function(x,P1=P1,P2=P2,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {
if (RS[n] == 0){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P2[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P2[n], log = TRUE) + dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}else{

Lp<-increment_log_prob(dbinom(0, size=1, prob=P2[n], log = TRUE) + dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }

###########################################################################################################
#Exposure Truncated Negative Binomial, Hurdle RS Negative Binomial
####################################################################### Define Custom Probability Function
dthnb<-function(x,P1=P1,P2=P2,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {

Lp<-increment_log_prob(dbinom(ifelse(RS[n]==0,1,0), size=1, prob=P2[n], log = TRUE));

if (RS[n] > 0){
Lp<- increment_log_prob(dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE)  -(1 - log(pnbinom(1, size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = FALSE)))  )  ;
}

}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }


###########################################################################################################
#Exposure Truncated Negative Binomial, RS Negative Binomial
####################################################################### Define Custom Probability Function
dtnb<-function(x,P1=P1,B1=B1,B2=B2,Mu1=Mu1,Mu2=Mu2, Exposure=Exposure,MaxExposure=MaxExposure, log=FALSE){
RS<-x

############################# Define log prob compilier
increment_log_prob <-function(ggg){
Lp<-ggg+Lp
}

################################################################################ Inflation model of Exposure
# Set Lp to Zero
N<-length(RS)
Lp<-0
for (n in 1:N) {
if(Exposure[n]==MaxExposure){
Lp<-increment_log_prob(dbinom(1, size=1, prob=P1[n], log = TRUE));

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )) ,log = TRUE)));

}else {

Lp<-increment_log_prob(dbinom(0, size=1, prob=P1[n], log = TRUE) + (dnbinom(Exposure[n], size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)-pnbinom((MaxExposure+1), size=(Mu1[n]/(B1[n])), prob=(1 / ( 1 + B1[n] )), log = TRUE)));

}}

################################################################################ ZINB Model of RS
for (n in 1:N) {
Lp<-increment_log_prob(dnbinom(RS[n], size=(Mu2[n]/(B2[n])), prob=(1 / ( 1 + B2[n] )), log = TRUE));

}
 if ( log==FALSE ) Lp <- exp(Lp)
return(Lp)
   }

