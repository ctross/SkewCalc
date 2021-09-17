model_code <-'
functions{
 real Mraw_index(vector r, vector t){
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
}

data{
 int N;
 int r[N];
 vector[N] t; 
 vector[N] t0;  
}

transformed data{
 int R;
 
 R = sum(r);
}

parameters{
 simplex[N] alpha;
 real<lower=0> gamma;
 real Concentration;
}

model{ 
 real T;
 real T_star;
 
 vector[N] t_hat;
 vector[N] t_hat_star;
 
 gamma ~ normal(1,0.01);
 Concentration ~ normal(0,2.5);
 
 T = sum(t-t0);
 t_hat = (t-t0)/T;
 
 T_star = sum(pow2(t,gamma) - pow2(t0,gamma));
 t_hat_star = (pow2(t,gamma) - pow2(t0,gamma))/T_star;
 
 alpha ~ dirichlet(t_hat_star*exp(Concentration));
 
 //r ~ multinomial(t_hat); 
 //r ~ multinomial(t_hat_star); 
 r ~ multinomial(alpha); 
}

generated quantities{
 real M_raw;
 real M;
 real M_raw_age;
 real M_age;
 
  { 
    real Bias;
    real T;
    real T_star;
 
    vector[N] t_hat;
    vector[N] t_hat_star;

    T = sum(t-t0);
    t_hat = (t-t0)/T;
 
    T_star = sum(pow2(t,gamma) - pow2(t0,gamma));
    t_hat_star = (pow2(t,gamma) - pow2(t0,gamma))/T_star;

    M_raw =  Mraw_index(to_vector(multinomial_rng(alpha,R)),t_hat);
    M_raw_age =  Mraw_index(to_vector(multinomial_rng(alpha,R)),t_hat_star);
 
    Bias = Mraw_index(to_vector(multinomial_rng(t_hat,R)),t_hat);
    M = M_raw - Bias;
 
    Bias = Mraw_index(to_vector(multinomial_rng(t_hat_star,R)),t_hat_star);
    M_age = M_raw_age - Bias;
 }
} 

'

#
