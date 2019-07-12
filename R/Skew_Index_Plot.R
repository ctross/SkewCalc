#' Simple posterior density plot for Stan results. 95 percent CI of posterior preditionsc in shaded orange, posterior median in white, and point estimate from sample in dark orange.
#'
#' @param Index Should plot be of "M" or of "Mraw"?
#' @param t M or Mraw should be adjusted for diminishing returns to age, TRUE or FALSE?
#' @param Save If desired, the plot can be save. Just set Save="FigueNameDesired", and R will export the plot as "FigureNameDesired.pdf".

#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index_stan(RS, Age)
#' skew_index_plot(Index="M", Age=FALSE, Save=FALSE)

skew_index_plot<-function(Index="M", Age=FALSE, SkewResults=StanReults, Save=FALSE){  
   
  if(Index=="M" && Age==FALSE){
  M <- extract(StanResults, pars="M")$M
  Point <- M_index(model_dat$r,model_dat$t) 
  }
  
  if(Index=="Mraw" && Age==FALSE){
  M <- extract(StanResults, pars="M_raw")$M_raw
  Point <- Mraw_index(model_dat$r,model_dat$t) 
  }
  
  if(Index=="M" && Age==TRUE){
  M <- extract(StanResults, pars="M_age")$M_age
  Point <- M_index_age(model_dat$r,model_dat$t,model_dat$t0) 
  }
  
  if(Index=="Mraw" && Age==TRUE){
  M <- extract(StanResults, pars="M_raw_age")$M_raw_age
  Point <- Mraw_index_age(model_dat$r,model_dat$t,model_dat$t0) 
  }
  
  q5  <- quantile(M, .05) 
  q95 <- quantile(M, .95)

  medx  <- median(M)
  x.dens  <- density(M)
  df.dens <- data.frame(x=x.dens$x, y=x.dens$y)
  
  dfx <- as.data.frame(x=M)
  dfP<<-dfx

p <- ggplot(data=dfx) + 
     geom_density(aes(x=M, y = ..density..), color = 'darkorange2') +
     geom_area(data = subset(df.dens, x >= q5 & x <= q95), 
              aes(x=x,y=y), fill='darkorange2', alpha=0.6) +
     geom_vline(xintercept=medx, color='#FFFFFF',size=2) +
     geom_vline(xintercept=Point, col="darkorange4")
 
print(p)
if(Save != FALSE){
 ggsave(file=paste0(Save,".pdf"), p,width=4,height=4) 
               }
            }
  

 
