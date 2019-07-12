model_code <- "
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
 
 Mraw = sqrt(N * sum(s)); 
 return Mraw;            
  }
 }

data{
 int N;
 int r[N];
 vector[N] t;	
 vector[N] t0;	
}

parameters{
 simplex[N] alpha;
 real<lower=0> gamma;
}

model{ 
 real T;
 real T_star;
 
 vector[N] t_hat;
 vector[N] t_hat_star;
 
 gamma ~ normal(0,1);
 
 T = sum((t-t0));
 t_hat = (t-t0)/T;
 
 T_star = sum((t^gamma - t0^gamma));
 t_hat_star = (t^gamma - t0^gamma)/T_star;
 
 alpha ~ dirichlet(rep_vector(1,N));
 
 r ~ multinomial(t_hat); 
 r ~ multinomial(t_hat_star); 
 r ~ multinomial(alpha); 
}

generated quantities{
 real M_raw;
 real M;
 real M_raw_age;
 real M_age;
 
{ real Bias;

 M_raw =  Mraw(alpha,t_hat);
 M_raw_age =  Mraw(alpha,t_hat_star);
 
 Bias = Mraw(to_vector(multinomial_rng((t-t0)/sum(t-t0),sum(r))),(t-t0));
 M = M_raw^2 - Bias^2;
 
 Bias = Mraw(to_vector(multinomial_rng((t^gamma - t0^gamma)/sum(t^gamma - t0^gamma),sum(r))),(t^gamma - t0^gamma));
 M_age = M_raw_age^2 - Bias^2;
}
}
"
