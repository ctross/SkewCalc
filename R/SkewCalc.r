# Code and Models By Cody T. Ross, SFI, UC Davis
 
####################################################################################### Run Model            
 SkewCalc <- function(RS,Exposure,
MaxExposure=max(Exposure),
Samples=1000,
MethodFit="ETNBShiftingShape",
MethodCI="HPDI",
P=0.95,
Prior="Default",
Updates="TRUE",
PlotResults="TRUE"
){

#############  Priors
if(Prior=="Default"){  
if(MethodFit=="ETZINBShiftingShape" | MethodFit=="ETZINBShiftingShapeScale"){
Priors<<-list( S1=c( 3.0,1.5),   # Scale Param, Log-scale, Constant, Age-Dist Model
               S2=c(-1.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
               S3=c(-1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
               Z1=c( 0.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
               Z2=c(-2.0,1.5),   # Zero inflation Param, Logit-scale, Constant or intecept, RS Model
               Z3=c( 0.0,1.5),   # Zero inflation Param, Logit-scale, slope on exposure, RS Model
               M1=c( 4.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
               M2=c(-2.0,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
               M3=c( 2.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )} else {
if(MethodFit=="ETHNBShiftingShape" | MethodFit=="ETHNBShiftingShapeScale"){
Priors<<-list( S1=c( 3.0,1.5),  # Scale Param, Log-scale, Constant, Age-Dist Model
             S2=c(-1.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
             S3=c( 1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
             Z1=c( 0.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
             Z2=c(-1.0,1.5),   # Zero inflation Param, Logit-scale, Constant or intecept, RS Model
             Z3=c( 0.0,1.5),   # Zero inflation Param, Logit-scale, slope on exposure, RS Model
             M1=c( 5.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
             M2=c( 1.0,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
             M3=c(-1.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )} 
              else {
if(MethodFit=="ETNBShiftingShape" | MethodFit=="ETNBShiftingShapeScale"){
Priors<<-list( S1=c( 4.0,1.5),   # Scale Param, Log-scale, Constant, Age-Dist Model
              S2=c(-2.0,1.5),   # Scale Param, Log-scale, Constant or Intercept, RS model
              S3=c(-1.0,1.5),   # Scale Param, Log-scale, Slope on Exposure, RS model
              Z1=c(-1.0,1.5),   # MaxExposure inflation Param, Logit-scale, Constant, Age-Dist Model
              M1=c( 5.0,1.5),   # Mu param, Log-scale, Constant, Age-Dist Model
              M2=c( 2.5,1.5),   # Mu param, Log-scale, Constant or Intercept, RS Model
              M3=c(-1.0,5.0)    # Mu param, Log-scale, Slope on Exposure, RS Model
             )}
             }} }
             else
             {Priors<<-Prior}
########################################################################################################
################################################################################# Define Models
flistShape <- alist(
RS ~ dzitnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

flistShapeScale <- alist(
RS ~ dzitnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

hlistShape <- alist(
RS ~ dthnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

hlistShapeScale <- alist(
RS ~ dthnb( P1,P2,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),
P2 <- logistic(Psi2 + Psi3*(Exposure)),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),
Psi2 ~ dnorm(Priors$Z2[1],Priors$Z2[2]),
Psi3 ~ dnorm(Priors$Z3[1],Priors$Z3[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

glistShape <- alist(
RS ~ dtnb( P1,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2*rep(1,length(Exposure))),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)

glistShapeScale <- alist(
RS ~ dtnb( P1,B1,B2,Mu1,Mu2, Exposure,MaxExposure) ,

P1 <- logistic(Psi1*rep(1,length(RS))),

B1 <- exp(Kappa1*rep(1,length(Exposure))),
B2 <- exp(Kappa2 + Kappa3*log(1+Exposure)),

Mu1 <-exp(Theta1*rep(1,length(Exposure))),
Mu2 <-exp(Theta2+ logistic(Theta3)*log(1+Exposure)),

Kappa1 ~ dnorm(Priors$S1[1],Priors$S1[2]),
Kappa2 ~ dnorm(Priors$S2[1],Priors$S2[2]),
Kappa3 ~ dnorm(Priors$S3[1],Priors$S3[2]),

Psi1 ~ dnorm(Priors$Z1[1],Priors$Z1[2]),

Theta1 ~ dnorm(Priors$M1[1],Priors$M1[2]),
Theta2 ~ dnorm(Priors$M2[1],Priors$M2[2]),
Theta3 ~ dnorm(Priors$M3[1],Priors$M3[2])
)


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# END LIBRARY

############################### Begin Function

if(Updates=="TRUE"){
print('Model has been built.')
print('Model is running, please be patient. The parameter estimation process may take several minutes with large datasets.')
}
############################################ Fit Model with MAP
if(MethodFit=="ETZINBShiftingShape"){
fitSkew <- map( flistShape , data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
 } else{
if(MethodFit=="ETZINBShiftingShapeScale"){
 fitSkew <- map(flistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
 }  else {
if(MethodFit=="ETHNBShiftingShape"){
fitSkew <- map( hlistShape , data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETHNBShiftingShapeScale"){
 fitSkew <- map(hlistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Psi2=Priors$Z2[1],Psi3=Priors$Z3[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETNBShiftingShape"){
 fitSkew <- map(glistShape, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Psi1=Priors$Z1[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else {
if(MethodFit=="ETNBShiftingShapeScale"){
 fitSkew <- map(glistShapeScale, data=list(RS=RS,Exposure=Exposure,MaxExposure=MaxExposure), start=list(Kappa1=Priors$S1[1],Kappa2=Priors$S2[1],Kappa3=Priors$S3[1],Psi1=Priors$Z1[1],Theta1=Priors$M1[1],Theta1=Priors$M2[1],Theta3=Priors$M3[1]) )
}
else(print('Error: Please Select a Valid Method for MethodFit')) }}} }}

if(Updates=="TRUE"){
print('Parameter estimation complete.')
print('Begining sampling from the posterior distribution.')
}

#################################################################################### Now Simulate Predictions
##########################################################################
######################################## Skew Calculation
M_1<-c()
Mc_1<-c()
GGG<-sample.qa.posterior(fitSkew,n=Samples)
FF<-Samples/10
Prog<-round(FF*c(1:10),0)

if(MethodFit=="ETZINBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(ZINBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETZINBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(ZINBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETHNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(HNBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETHNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(HNBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(NBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
 else{
if(MethodFit=="ETNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx<-M_Mc(NBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS)))
M_1[i] <-xxxx[1]
Mc_1[i] <-xxxx[2]
} }
else(print('Error: Please Select a Valid Method for MethodFit')) }}}}}

if(MethodCI=="HPDI"){
  RES<-matrix(NA,ncol=4,nrow=4)

  RES[1,] <-c(M_Mc(c(RS,Exposure))[1],NA,NA,NA)
  RES[2,] <-c(mean(M_1),median( M_1),HPDI(M_1,P)[1],HPDI(M_1,P)[2])

  RES[3,] <-c(M_Mc(c(RS,Exposure))[2],NA,NA,NA)
  RES[4,] <-c(mean(Mc_1),median( Mc_1),HPDI(Mc_1,P)[1],HPDI(Mc_1,P)[2])

  colnames(RES)<-c("MAP","Median","HPDI-L","HPDI-H")
  rownames(RES)<-c("Sample M","Population M","Sample Mc","Population Mc" )

 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning(paste('The Sample Prediction Of M or Mc Is Not Contained in the Models Central' ,(P*100) ,'Percent Posterior Predictive Density.'))
 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning('The Selected Exposure-Truncated Negative-Binomial Model Might Not Be a Good Approximating Model For This Sample.')

} else{
if(MethodCI=="PCI"){

  RES<-matrix(NA,ncol=4,nrow=4)

  RES[1,] <-c(M_Mc(c(RS,Exposure))[1],NA,NA,NA)
  RES[2,] <-c(mean(M_1),median( M_1),PCI(M_1,P)[1],PCI(M_1,P)[2])

  RES[3,] <-c(M_Mc(c(RS,Exposure))[2],NA,NA,NA)
  RES[4,] <-c(mean(Mc_1),median( Mc_1),PCI(Mc_1,P)[1],PCI(Mc_1,P)[2])

  colnames(RES)<-c("MAP","Median","PCI-L","PCI-H")
  rownames(RES)<-c("Sample M","Population M","Sample Mc","Population Mc" )

 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning(paste('The Sample Prediction Of M or Mc Is Not Contained in the Models Central' ,(P*100) ,'Percent Posterior Predictive Density.'))
 if((RES[1,1]<RES[2,3] | RES[1,1]>RES[2,4]) |(RES[3,1]<RES[4,3] | RES[3,1]>RES[4,4])) warning('The Selected The Exposure-Truncated Zero-Inflated Negative-Binomial Model Might Not Be a Good Approximating Model For This Sample.')

 }
else ( print('Error in Specifying Confidence Interval Method'))
 }
 if(PlotResults=="TRUE"){
print('Plotting requires extra computation, please be patient.')
######################################## Sample Data
library(RColorBrewer)
Samples<-100
FF<-Samples/10
Prog<-round(FF*c(1:10),0)
plotRS<-c()
plotExposure<-c()

if(MethodFit=="ETZINBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-ZINBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
} }

 else{
if(MethodFit=="ETZINBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-ZINBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETHNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-HNBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETHNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-HNBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETNBShiftingShape"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-NBShape_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }
 else{
if(MethodFit=="ETNBShiftingShapeScale"){
for( i in 1:Samples){
if(Updates=="TRUE"){
for(j in 1:10){
if(i==Prog[j])print(paste('Sampling is',(Prog[j]/Samples)*100, 'percent complete.') )
}  }

xxxx2<-NBShapeScale_rng(i=1, G=GGG, MaxExposure, N=length(RS))
plotRS <-c(plotRS,xxxx2[1:length(RS)])
plotExposure <-c(plotExposure,xxxx2[(1+length(RS)):(length(RS)*2)])
}  }

else(print('Error: Please Select a Valid Method for MethodFit'))  }}} } }

 
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

 ######################################################################## Plot it
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths =c(2,2) , heights=c(1,1))

   jet.colors <-colorRampPalette(c("white",brewer.pal(9,"YlOrRd")))
 MaxRS<-max(c(max(RS,na.rm=T),max(plotRS,na.rm=T)),na.rm=T)
 smoothScatter(RS~Exposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                 colramp = jet.colors,transformation = function(x) x^.45,
                 nbin=200,  xlab="Sample Age Distribution",
                 ylab="Sample RS Distribution")

 smoothScatter(plotRS~plotExposure,xlim=c(0,MaxExposure),ylim=c(0,MaxRS),
                colramp = jet.colors,transformation = function(x) x^.45,
                nbin=200, xlab="Population Age Distribution (Predictions)",
                ylab="Predicted Posterior Population RS Distribution")

 }


if(MethodCI=="HPDI"){
print(paste('Confidence Intervals are: ',P*100,'percent Highest Posterior Density Intervals.'))
} else{
if(MethodCI=="PCI"){
print(paste('Confidence Intervals are: Central', P*100,'percent Posterior Credibility Intervals.'))
 }
else ( print('Error in Specifying Confidence Interval Method'))
 }
 print(precis(fitSkew))
 print(RES)
  if(PlotResults=="TRUE"){
return(list(RES=RES,Precis=precis(fitSkew),P_RS=plotRS,P_Exposure=plotExposure))
} else{
return(list(RES=RES,Precis=precis(fitSkew)))
}
 }
 

################################################################################ Example with Fake Data
################ Set Parameters
#N<-100
#MaxExposure<-45

#Kappa1 <- 3
#Kappa2 <- -2
#Psi1 <- -1
#Theta1 <- 5
#Theta2 <- 1.5
#Theta3 <- .5

#scrapExposure<-c()
#scrapRS<-c()

################################################ Split Exposure RNG into point process at MaxExposure and the NB
#Nscrap<-round(logistic(Psi1)*N);
#Npp <- 1;
#while (Npp < Nscrap) {
#Npp <- Npp + 1;
#}
################################################################################ Model of Exposure
#  for( n in 1:Npp){
#       scrapExposure[n]<-MaxExposure;
#       }
#  for( n in (Npp+1):N){
#      Ticker <- 1;
#      while (Ticker == 1) {
#        Nscrap <-rnbinom(1, size=(exp(Theta1)/(exp(Kappa1))), prob=(1 / ( 1 + (exp(Kappa1)) )))
#        if(Nscrap<=MaxExposure){
#        Ticker<-0;}
#        else{
#        Ticker<-1;}}

#       scrapExposure[n]<-Nscrap;
#       }
################################################################################ Model of RS
#   for( n in 1:N){
#   scrapRS[n]<- rnbinom( 1,size=(exp(Theta2+ logistic(Theta3)*log(1+scrapExposure[n]))/ (exp(Kappa2))),prob=(1 / ( 1 + (exp(Kappa2)) )));
#         }
################################################################################ Run SkewCalc
# RES<-SkewCalc(scrapRS,scrapExposure, MethodFit="ETNBShiftingShape")
# Parameters should be recovered
 
################################################################################ Example with Real Data
# data(Kipsigis)
# d<-Kipsigis
# RS<-round(d$soff,0)
# Exposure<-d$age
# Exposure<-round(Exposure-15,0)
# Exposure<-ifelse(Exposure>45,45,Exposure)
# MaxExposure<-max(Exposure)

# RES<-SkewCalc(RS,Exposure,Samples=100, MethodFit="ETNBShiftingShape")
    
 
 
 
