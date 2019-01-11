model_code <-"
functions{
real M(int[] r, int[] t1, int[] t0, int N){
 real R; 
 vector[N] tt; 
 real T; 
 vector[N] si; 
 real S;  
 real M1;                      

 R = sum(r);
 
 for(i in 1:N) 
 tt[i] = t1[i] - t0[i];
 
 T = sum(tt);  
 
 for(i in 1:N)
 si[i] = ((r[i]/R)-(tt[i]/T))^2;  
 
 S = sum(si);  
 
 M1 = sqrt(N * S)*sqrt(R/N);  
 return M1;          
 }
 
 real Mraw(int[] r, int[] t1, int[] t0, int N){
 real R; 
 vector[N] tt; 
 real T; 
 vector[N] si; 
 real S; 
 real Mraw1;                       

 R = sum(r);
 
 for(i in 1:N) 
 tt[i] = t1[i] - t0[i];
 
 T = sum(tt);  
 
 for(i in 1:N)
 si[i] = ((r[i]/R)-(tt[i]/T))^2;  
 
 S = sum(si);  
 
 Mraw1 = sqrt(N * S); 
 return Mraw1;            
 }
 
real M_age(int[] r, int[] t1, int[] t0, real Beta, int N){
 real R; 
 vector[N] tt; 
 real T; 
 vector[N] si; 
 real S;  
 real M1;                      

 R = sum(r);
 
 for(i in 1:N) 
 tt[i] = t1[i]^Beta - t0[i]^Beta;
 
 T = sum(tt);  
 
 for(i in 1:N)
 si[i] = ((r[i]/R)-(tt[i]/T))^2;  
 
 S = sum(si);  
 
 M1 = sqrt(N * S)*sqrt(R/N);  
 return M1;          
 }
 
 real Mraw_age(int[] r, int[] t1, int[] t0, real Beta, int N){
 real R; 
 vector[N] tt; 
 real T; 
 vector[N] si; 
 real S; 
 real Mraw1;                       

 R = sum(r);
 
 for(i in 1:N) 
 tt[i] = t1[i]^Beta - t0[i]^Beta;
 
 T = sum(tt);  
 
 for(i in 1:N)
 si[i] = ((r[i]/R)-(tt[i]/T))^2;  
 
 S = sum(si);  
 
 Mraw1 = sqrt(N * S); 
 return Mraw1;            
 }


}

data{
 int N;
 int T0[N];
 int T1[N];
 int RS[N];
}

transformed data{
 int Z; 
 int Q;
 int Fail;
 
 Fail = 5;    
 
 Z = max(T1);
 Q = min(T1);
}

parameters{
 real<lower=0, upper=1> Mu_T0;
 real<lower=0> B_T0;
 
 real<lower=0> Mu_T1;
 real<lower=0> B_T1;
 
 real<lower=0> A_RS;  
 real<lower=0, upper=1> Beta;
 real<lower=0> B_RS;
}

model{
 vector[N] Mu_RS;
 
//# Model Exposure Duration
 Mu_T0 ~ beta(4, 1);
 B_T0 ~ normal(0,10);

 for(i in 1:N)
 T0[i] ~ beta_binomial(T1[i],Mu_T0*(B_T0+1),(1-Mu_T0)*(B_T0+1))T[0,T1[i]-1];
 
//# Model Age of Death or Censor 
 Mu_T1 ~ gamma(1, 0.05);
 B_T1 ~ normal(0,10);

 for(i in 1:N)
 T1[i] ~ neg_binomial(Mu_T1*B_T1,B_T1)T[Q,Z];
 
//# Model RS 
 A_RS ~ normal(0,10);
 B_RS ~ normal(0,10);
 
 Beta ~ beta(2,2);
 
 for(i in 1:N){
   Mu_RS[i] = A_RS*(T1[i]^Beta - T0[i]^Beta);
               }
               
 RS ~ neg_binomial(Mu_RS*B_RS,B_RS);
}

generated quantities{
 int P_T1[N];
 int P_T0[N];
 int P_RS[N];
 
 vector[4] Mraw_index;
 vector[4] M_index; 
 {
 int y;
 
 for(i in 1:N){
    y = 1;    
    while(y < Fail){
     P_T1[i] = neg_binomial_rng(Mu_T1*B_T1,B_T1);
     if( P_T1[i] <= Z && P_T1[i]>=Q )
      {y = Fail;} else {
      y = y+1;
      P_T1[i] = Z;}
    }

   y = 1;    
   while(y < Fail){      
    P_T0[i] = beta_binomial_rng(P_T1[i],Mu_T0*(B_T0+1),(1-Mu_T0)*(B_T0+1));
    if( P_T0[i] < P_T1[i] )
      {y = Fail;} else {
      y = y+1;
      P_T0[i] = 0;}
    }
    
    P_RS[i] = neg_binomial_rng((A_RS*(P_T1[i]^Beta - P_T0[i]^Beta))*B_RS,B_RS);    
              }  
 }             
 
 Mraw_index[1] = Mraw(RS,T1,T0,N);
 M_index[1] = M(RS,T1,T0,N);
 
 Mraw_index[2] = Mraw_age(RS,T1,T0,Beta,N);
 M_index[2] = M_age(RS,T1,T0,Beta,N);
 
 Mraw_index[3] = Mraw(P_RS,P_T1,P_T0,N);
 M_index[3] = M(P_RS,P_T1,P_T0,N);
              
 Mraw_index[4] = Mraw_age(P_RS,P_T1,P_T0,Beta,N);
 M_index[4] = M_age(P_RS,P_T1,P_T0,Beta,N); 
}
"
