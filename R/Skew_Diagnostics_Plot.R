#' Check the fit of the predictive model to the sample data
#'
#' @param r A vector of RS values.
#' @param t A vector of ages at death or last census.
#' @param t A vector of ages at first census. Defaults to zero.
#' @param Save If desired, the plot can be save. Just set Save="FigueNameDesired", and R will export the plot as "FigureNameDesired.pdf".
#' @examples
#' set.seed(1)
#' RS = rpois(100, 5) 
#' Age = rpois(100, 45)
#' M_index_stan(RS, Age)
#' skew_index_plot(Index="M", Age=FALSE, Save=FALSE)

skew_diagnostics_plot<-function(RS,T1,T0=0,SkewResults=StanReults, save=FALSE){  
  P_RS<-c(extract(StanResults, pars="P_RS")$P_RS)
  P_T1<-extract(StanResults, pars="P_T1")$P_T1
  P_T0<-extract(StanResults, pars="P_T0")$P_T0
  P_A <- c(P_T1 - P_T0)
  A <- T1-T0

  SP_A <- sample(P_A,length(A)*100)
  SP_RS <- sample(P_RS,length(RS)*100)

  df1 <- data.frame(RS=c(RS,P_RS), Age=c(A,P_A),ID=c(rep("Sample",length(RS)),rep("Simulated",length(P_RS))))
  df2 <- df1[which(df1$ID=="Sample"),] 
 
  grob1 <- ggplot(df1, aes(RS, fill = ID, colour = ID)) +
       geom_density(alpha = 0.3,adjust = 2,size=1) + 
       scale_colour_manual(values = c("black", "darkorange2")) + 
       scale_fill_manual(values = c("black", "darkorange2")) 
  
  grob2 <- ggplot(df1, aes(Age, fill = ID, colour = ID)) +
       geom_density(alpha = 0.3,adjust = 2,size=1) + 
       scale_colour_manual(values = c("black", "darkorange2")) + 
       scale_fill_manual(values = c("black", "darkorange2")) 
  
   
  grob3 <- ggplot(df1, aes(x=Age, y=RS) )  +
       geom_hex(bins = 25) +
       scale_fill_gradient(low = "gainsboro", high = "darkorange2") + 
       geom_point(data = df2)
  
grid.arrange(grob1, grob2, grob3, nrow = 1)  

if(save != FALSE){
 g <- arrangeGrob(grob1, grob2, grob3, nrow = 1) 
 ggsave(file=paste0(save,".pdf"), g,width=12,height=4) 
               }
            }
  






















