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
  M <- extract(StanResults, pars="M_index")$M_index
  Point <- median(M[,1])
  M <- M[,3]
  }
  
  if(Index=="Mraw" && Age==FALSE){
  M <- extract(StanResults, pars="Mraw_index")$Mraw_index
  Point <- median(M[,1])
  M <- M[,3]
  }
  
  if(Index=="M" && Age==TRUE){
  M <- extract(StanResults, pars="M_index")$M_index
  Point <- median(M[,2])
  M <- M[,4]
  }
  
  if(Index=="Mraw" && Age==TRUE){
  M <- extract(StanResults, pars="Mraw_index")$Mraw_index
  Point <- median(M[,2])
  M <- M[,4]
  }
  
  q2.5  <- quantile(M, .025) 
  q97.5 <- quantile(M, .975)

  medx  <- median(M)
  x.dens  <- density(M)
  df.dens <- data.frame(x=x.dens$x, y=x.dens$y)
  
  dfx <- as.data.frame(x=M)

p <- ggplot(data=dfx) + 
     geom_density(aes(x=M, y = ..density..), color = 'darkorange2') +
     geom_area(data = subset(df.dens, x >= q2.5 & x <= q97.5), 
              aes(x=x,y=y), fill='darkorange2', alpha=0.6) +
    geom_vline(xintercept=Point, col="darkorange4") +
    geom_vline(xintercept=medx, color='#FFFFFF',size=2)
 
print(p)
if(Save != FALSE){
 ggsave(file=paste0(Save,".pdf"), p,width=4,height=4) 
               }
            }
  

 