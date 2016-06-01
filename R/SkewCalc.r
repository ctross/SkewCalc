#
# If you don't have Stan, download and install R-Tools 
# https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows
#
# Then follow instructions https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started, or just run
# source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
# install_rstan()#
#
####################
 N25<- function(x){
 length(x[which(x>=25)])
 }

 RS0A25<- function(x,y){
 dat<-data.frame(x,y)
 dat2<-dat[which(dat[,1]>=25),]
 length(which(dat2[,2]==0))
 }

 PercRS0A25<- function(x,y){
 dat<-data.frame(x,y)
 dat2<-dat[which(dat[,1]>=25),]
 round(length(which(dat2[,2]==0))/length(dat2[,2]),2)
 }

 N45<- function(x){
 length(x[which(x>=45)])
 }

 RS0A45<- function(x,y){
 dat<-data.frame(x,y)
 dat2<-dat[which(dat[,1]>=45),]
 length(which(dat2[,2]==0))
 }

 PercRS0A45<- function(x,y){
 dat<-data.frame(x,y)
 dat2<-dat[which(dat[,1]>=45),]
 round(length(which(dat2[,2]==0))/length(dat2[,2]),2)
 }
##################################

N0<- function(x){
 length(x)
 }

RSmean<- function(y){
 mean(y)
 }

RSsum<- function(y){
 sum(y)
 }

Exposuresum<- function(x){
 sum(x)
 }

RSsd<- function(y){
 sd(y)
 }

RSvar<- function(y){
 var(y)
 }

RScv<- function(y){
 sd(y)/mean(y)
 }

RSofs<- function(y){
 var(y)/(mean(y)*mean(y))
 }

Gini<-function (y, corr = FALSE, na.rm = TRUE)
{
    if (!na.rm && any(is.na(y)))
        return(NA_real_)
    y <- as.numeric(na.omit(y))
    n <- length(y)
    y <- sort(y)
    G <- sum(y * 1L:n)
    G <- 2 * G/sum(y) - (n + 1L)
    if (corr)
        G/(n - 1L)
    else G/n
}

RateGini<-function (y,x, corr = FALSE, na.rm = TRUE)
{
y<-y/x
    if (!na.rm && any(is.na(y)))
        return(NA_real_)
    y <- as.numeric(na.omit(y))
    n <- length(y)
    y <- sort(y)
    G <- sum(y * 1L:n)
    G <- 2 * G/sum(y) - (n + 1L)
    if (corr)
        G/(n - 1L)
    else G/n
}

Bindex <- function(ki,ni) {   #ki=rs,ni=exposure
	Nt <- sum(ni)
	Nbar <- Nt / max(ni)
	K <- sum(ki)
	if(K>0){
		B = sum((ki / K - ni / Nt)^2) - (1 / K) * (1 - 1 / Nbar)
	} else {
		B = 0
	}
	B
}

Mraw_index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	C <- sqrt(N * S)
  C
}

M_index <- function(ki,ni) {
	N <- length(ki)
	K <- sum(ki)
	fi <- ni
	SumF <- sum(fi)
	si <- ((ki/K)-(fi/SumF))^2
	S <- sum(si)
	D1 <- sqrt(N * S)*sqrt(mean(ki))
  D1
}

 library(Cairo)

 MergeRes<-function(Exposure,RS,NHM=FALSE){
                            if(NHM==FALSE){
 Exposure<-ifelse(Exposure>60,60,Exposure)
 Exposure <- Exposure - 11
         }
   SkewCalc(RS,Exposure, Samples=1000, Warmup=500, Chains=1, Refresh=1, Code="GRF")

  c(
  N25(Exposure),RS0A25(Exposure,RS),PercRS0A25(Exposure,RS),N45(Exposure),RS0A45(Exposure,RS),PercRS0A45(Exposure,RS),
  N0(RS),RSsum(RS),Exposuresum(Exposure),RSmean(RS),RSsd(RS),RSvar(RS),RScv(RS), RSofs(RS),Gini(RS),
  RateGini(RS,Exposure),Bindex(RS,Exposure),Mraw_index(RS,Exposure),M_index(RS,Exposure),
  c(SkewResults$SkewFit),    c(SkewResults$SkewGammaHyperParams)
  )
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
 StanResults <- stan(model_code=skew_code_GRF, data=model_dat, thin=1, iter=Samples, warmup=Warmup, chains=Chains, refresh=Refresh,init=0)
 MMC<-extract(StanResults, pars="Mraw_M")$Mraw_M
 MMC<-MMC[which(MMC[,1] < 9999),]
 MMC<-MMC[which(MMC[,2] < 9999),]
 
SkewGammaHyperParams<-matrix(NA,nrow=2,ncol=2)
SkewGammaHyperParams[1,1]<-fitdistr(MMC[,1],"gamma")$estimate[1]
SkewGammaHyperParams[1,2]<-fitdistr(MMC[,1],"gamma")$estimate[2]
SkewGammaHyperParams[2,1]<-fitdistr(MMC[,2],"gamma")$estimate[1]
SkewGammaHyperParams[2,2]<-fitdistr(MMC[,2],"gamma")$estimate[2]

 colnames(SkewGammaHyperParams)<-c("Shape","Scale")
 rownames(SkewGammaHyperParams)<-c("Mraw","M")

 SkewFit<-matrix(NA,nrow=4,ncol=3)
 SkewFit[1,1]<-Mraw_index(RS,Exposure)
 SkewFit[2,1]<-M_index(RS,Exposure)
 SkewFit[3,1]<-median(MMC[,1])
 SkewFit[4,1]<-median(MMC[,2])
 SkewFit[3,2]<-HPDI(MMC[,1])[1]
 SkewFit[4,2]<-HPDI(MMC[,2])[1]
 SkewFit[3,3]<-HPDI(MMC[,1])[2]
 SkewFit[4,3]<-HPDI(MMC[,2])[2]

 colnames(SkewFit)<-c("Estimate","Lower_90%_HPDI","Upper_90%_HPDI")
 rownames(SkewFit)<-c("Mraw_Point_Estimate","M_Point_Estimate","Mraw_Posterior_Estimate","M_Posterior_Estimate")

 print(SkewFit)
 SkewResults<<-list(SkewFit=SkewFit,StanResults=StanResults,SkewGammaHyperParams=SkewGammaHyperParams)}
 else{   if(Code=="Fast"){
 StanResults <- stan(model_code=skew_code_Fast, data=model_dat, thin=1, iter=Samples, warmup=Warmup, chains=Chains, refresh=Refresh,init=0)
       MMC<-extract(StanResults, pars="Mraw_M")$Mraw_M
        MMC<-MMC[which(MMC[,1] < 9999),]
 MMC<-MMC[which(MMC[,2] < 9999),]
       
SkewGammaHyperParams<-matrix(NA,nrow=2,ncol=2)
SkewGammaHyperParams[1,1]<-fitdistr(MMC[,1],"gamma")$estimate[1]
SkewGammaHyperParams[1,2]<-fitdistr(MMC[,1],"gamma")$estimate[2]
SkewGammaHyperParams[2,1]<-fitdistr(MMC[,2],"gamma")$estimate[1]
SkewGammaHyperParams[2,2]<-fitdistr(MMC[,2],"gamma")$estimate[2]

 colnames(SkewGammaHyperParams)<-c("Shape","Scale")
 rownames(SkewGammaHyperParams)<-c("Mraw","M")

 SkewFit<-matrix(NA,nrow=4,ncol=3)
 SkewFit[1,1]<-Mraw_index(RS,Exposure)
 SkewFit[2,1]<-M_index(RS,Exposure)
 SkewFit[3,1]<-median(MMC[,1])
 SkewFit[4,1]<-median(MMC[,2])
 SkewFit[3,2]<-HPDI(MMC[,1])[1]
 SkewFit[4,2]<-HPDI(MMC[,2])[1]
 SkewFit[3,3]<-HPDI(MMC[,1])[2]
 SkewFit[4,3]<-HPDI(MMC[,2])[2]

 colnames(SkewFit)<-c("Estimate","Lower_90%_HPDI","Upper_90%_HPDI")
 rownames(SkewFit)<-c("Mraw_Point_Estimate","M_Point_Estimate","Mraw_Posterior_Estimate","M_Posterior_Estimate")

 print(SkewFit)
 SkewResults<<-list(SkewFit=SkewFit,StanResults=StanResults,SkewGammaHyperParams=SkewGammaHyperParams)}
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
  Exposure<-ifelse(Exposure>60,60,Exposure)
  Exposure <- Exposure - 11
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

########### Sample Data
SampleAge<-c(27,45,42,24,24,21,17,25,13,28,25,47,50,43,29,17,31,45,34,21,17,36,33,41,37,46,42,15,46,45,35,30,11,44,44,45,28,46,25,15,18,57,30,45,46,40,58,24,29,28,24,25,24,44,20,23,18,29,32,30,26,26,32)
SampleExposure<-ifelse(SampleAge>44,45,SampleAge)
SampleRS<-c(10,9,17,2,5,9,0,5,2,2,6,20,15,6,2,5,12,14,1,3,2,11,3,15,9,7,10,1,16,11,15,12,0,26,9,5,2,9,8,1,1,5,14,27,5,12,33,7,5,3,4,3,4,13,4,4,1,4,30,7,6,3,4)
