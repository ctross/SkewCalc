model_code_grf <-'
functions{
 real Mraw(vector r, vector t){
   int N = rows(r);

   real R; 
   real T; 
   vector[N] s; 
   real Mraw;                       

   R = sum(r);
   T = sum(t);  

   for(i in 1:N)
   s[i] = ((r[i]/R)-(t[i]/T))^2;  
 
   Mraw = N * sum(s); 
   return Mraw;            
  }

 vector pow2(vector x, real y){
   vector[rows(x)] z;
   for (i in 1:rows(x)) 
   z[i] = x[i]^y;
   return(z);
  }

 //### Cholesky Factor of Covariance Matrix from a Gaussian Process 
  matrix GP(int K, real C, real D, real S){
   matrix[K,K] Rho;                       
   real KR;                               
   KR = K;                             

   for(i in 1:(K-1)){
   for(j in (i+1):K){  
    Rho[i,j] = C * exp(-D * ( (j-i)^2 / KR^2) );           
    Rho[j,i] = Rho[i,j];                       
    }}

   for (i in 1:K){
    Rho[i,i] = 1;                               
    }

   return S*cholesky_decompose(Rho);
   }
}

data{
 int N;
 int A;
 int r[N];
 int t[N]; 
 int t0[N];  
}

transformed data{
 int R;
 vector[N] t_real; 
 vector[N] t0_real;  
 
 for( i in 1:N){
 t_real[i] = t[i];
 t0_real[i] = t0[i];
 }

 R = sum(r);
}

parameters{
 simplex[N] alpha;
 real Concentration;
 vector[A] theta_raw;
 real<lower=0> S;
 real<lower=0> D;
 real<lower=0, upper=1> C;
 real Mu;
}

model{ 
 real T;
 real T_star;

 int t0p1[N];

 vector[N] t_eff;
 vector[N] t_hat;
 vector[N] t_hat_star;
 vector[A] theta;
 
 Concentration ~ normal(0,2.5);

 //# Priors for Gaussian Process controlers
  S ~ exponential(1);
  D ~ exponential(1);
  C ~ beta(12, 2);
  Mu ~ normal(0, 1);
 
  theta_raw ~ normal(0,1);
  theta = Mu + GP(A, C, D, S)*theta_raw;
  
 for(i in 1:N){
  t_eff[i] = 0;
  t0p1[i] = t0[i] + 1;
   for(a in t0p1[i]:t[i]){
    t_eff[i] += inv_logit(theta[a]);
    }}
 
 T = sum(t_real-t0_real);
 t_hat = (t_real-t0_real)/T;
 
 T_star = sum(t_eff);
 t_hat_star = (t_eff)/T_star;
 
 alpha ~ dirichlet(t_hat*exp(Concentration));
 
 r ~ multinomial(t_hat); 
 r ~ multinomial(t_hat_star); 
 r ~ multinomial(alpha); 
}

generated quantities{
 real M_raw;
 real M;
 real M_raw_age;
 real M_age;
 vector[A] theta_normalized;

 theta_normalized = cumulative_sum( inv_logit(Mu + GP(A, C, D, S)*theta_raw ));

  { 
    real Bias; 
    real T;
    real T_star;
    vector[A] theta;
 
    int t0p1[N];
    vector[N] t_eff;
    vector[N] t_hat;
    vector[N] t_hat_star;

    theta = Mu + GP(A, C, D, S)*theta_raw;

    for(i in 1:N){
     t_eff[i] = 0;
     t0p1[i] = t0[i] + 1;
     for(a in t0p1[i]:t[i]){
     t_eff[i] += inv_logit(theta[a]);
      }}
 
    T = sum(t_real-t0_real);
    t_hat = (t_real-t0_real)/T;
 
    T_star = sum(t_eff);
    t_hat_star = (t_eff)/T_star;

    M_raw =  Mraw(to_vector(multinomial_rng(alpha,R)),t_hat);
    M_raw_age =  Mraw(to_vector(multinomial_rng(alpha,R)),t_hat_star);
 
    Bias = Mraw(to_vector(multinomial_rng(t_hat,R)),t_hat);
    M = M_raw - Bias;
 
    Bias = Mraw(to_vector(multinomial_rng(t_hat_star,R)),t_hat_star);
    M_age = M_raw_age - Bias;
 }
} 

'

#
