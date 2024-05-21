# JOINT SPECIFICATION

jointmodel <- "
// BI-EXPONENTIAL SPECIFICATION
functions {
  vector nonlinear_predictor(int[] ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  int n;
  int nbeta;
  vector[N_M] y_M;
  vector[N_F] y_F;
  int<lower=1,upper=n_M> ID_M1[N_M];
  int<lower=1,upper=n> ID_M2[N_M];
  int<lower=1,upper=n_F> ID_F1[N_F];
  int<lower=1,upper=n> ID_F2[N_F];
  vector[N_M] times_M;
  vector[N_F] times_F;
  matrix[n,nbeta] X;
  vector[n] Time;
  vector[n] status_1;
  vector[n] status_2;
}
parameters {
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  cov_matrix[3] Omega_M;
  cov_matrix[3] Omega_F;
  vector[nbeta] beta_1;
  vector[nbeta] beta_2;
  vector[2] alphaB_1;
  vector[2] alphaG_1;
  vector[2] alphaD_1;
  vector[2] alphaB_2;
  vector[2] alphaG_2;
  vector[2] alphaD_2;
  real gamma_1;
  real gamma_2;
  real<lower=0> phi_1;
  real<lower=0> phi_2;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
model {
  vector[n] loghaz_1;
  vector[n] loghaz_2;
  vector[n] cumHaz_1;
  vector[n] cumHaz_2;

  // BI-EXPONENTIAL SPECIFICATION
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y_M | nonlinpred_M, sqrt(sigma2_M));
  target += normal_lpdf(y_F | nonlinpred_F, sqrt(sigma2_F));

  // COMPETING RISKS LOG-LIKELIHOOD
  for(i in 1:n){
       // Hazard functions
       loghaz_1[i] = log(phi_1) + (phi_1-1)*log(Time[i]) + gamma_1 + X[i,]*beta_1 + 
                       lBi[i,]*alphaB_1 + lGi[i,]*alphaG_1 + lDi[i,]*alphaD_1;
       loghaz_2[i] = log(phi_2) + (phi_2-1)*log(Time[i]) + gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2;

       // Cumulative hazard functions
       cumHaz_1[i] = pow(Time[i], phi_1) * exp( gamma_1 + X[i,]*beta_1 + 
                       lBi[i,]*alphaB_1 + lGi[i,]*alphaG_1 + lDi[i,]*alphaD_1 );
       cumHaz_2[i] = pow(Time[i], phi_2) * exp( gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2 );
       
       target += status_1[i]*loghaz_1[i] + status_2[i]*loghaz_2[i] - cumHaz_1[i] - cumHaz_2[i];
  }

  // LOG-PRIORS
  // Longitudinal fixed effects
  target += normal_lpdf(theta_M | 0, 10);
  target += normal_lpdf(theta_F | 0, 10);
   
  // Residual error variance
  target += cauchy_lpdf(sigma2_M | 0, 5);
  target += cauchy_lpdf(sigma2_F | 0, 5);
 
  // Random-effects variance-covariance matrices
  target += inv_wishart_lpdf(Omega_M | 4, diag_matrix(rep_vector(1,3)));
  target += inv_wishart_lpdf(Omega_F | 4, diag_matrix(rep_vector(1,3)));
  
  // Random-effects
  for(i in 1:n_M){ target += multi_normal_lpdf(bi_M[i,1:3] | rep_vector(0,3), Omega_M); }
  for(i in 1:n_F){ target += multi_normal_lpdf(bi_F[i,1:3] | rep_vector(0,3), Omega_F); }
  
  // Survival fixed effects
  target += normal_lpdf(beta_1 | 0, 10);  
  target += normal_lpdf(beta_2 | 0, 10);
  target += normal_lpdf(gamma_1 | 0, 10);  
  target += normal_lpdf(gamma_2 | 0, 10);
 
  // Association parameters
  target += normal_lpdf(alphaB_1 | 0, 10);
  target += normal_lpdf(alphaB_2 | 0, 10);
  target += normal_lpdf(alphaG_1 | 0, 10);
  target += normal_lpdf(alphaG_2 | 0, 10);
  target += normal_lpdf(alphaD_1 | 0, 10);
  target += normal_lpdf(alphaD_2 | 0, 10);

  // Shape parameters (Weibull hazard)
  target += cauchy_lpdf(phi_1 | 0, 1);
  target += cauchy_lpdf(phi_2 | 0, 1);
  
}

"


jointmodelsurv <- "
// BI-EXPONENTIAL SPECIFICATION
functions {
  vector nonlinear_predictor(int[] ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  int n;
  int nbeta;
  vector[N_M] y_M;
  vector[N_F] y_F;
  int<lower=1,upper=n_M> ID_M1[N_M];
  int<lower=1,upper=n> ID_M2[N_M];
  int<lower=1,upper=n_F> ID_F1[N_F];
  int<lower=1,upper=n> ID_F2[N_F];
  vector[N_M] times_M;
  vector[N_F] times_F;
  matrix[n,nbeta] X;
  vector[n] Time;
  vector[n] status_2;
}
parameters {
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  cov_matrix[3] Omega_M;
  cov_matrix[3] Omega_F;
  vector[nbeta] beta_2;
  vector[2] alphaB_2;
  vector[2] alphaG_2;
  vector[2] alphaD_2;
  real gamma_2;
  real<lower=0> phi_2;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
model {
  vector[n] loghaz_2;
  vector[n] cumHaz_2;

  // BI-EXPONENTIAL SPECIFICATION
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y_M | nonlinpred_M, sqrt(sigma2_M));
  target += normal_lpdf(y_F | nonlinpred_F, sqrt(sigma2_F));

  // SURVIVAL LOG-LIKELIHOOD
  for(i in 1:n){
       // Hazard functions
       loghaz_2[i] = log(phi_2) + (phi_2-1)*log(Time[i]) + gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2;

       // Cumulative hazard functions
       cumHaz_2[i] = pow(Time[i], phi_2) * exp( gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2 );
       
       target += status_2[i]*loghaz_2[i] - cumHaz_2[i];
  }

  // LOG-PRIORS
  // Longitudinal fixed effects
  target += normal_lpdf(theta_M | 0, 10);
  target += normal_lpdf(theta_F | 0, 10);
   
  // Residual error variance
  target += cauchy_lpdf(sigma2_M | 0, 5);
  target += cauchy_lpdf(sigma2_F | 0, 5);

  // Random-effects variance-covariance matrices
  target += inv_wishart_lpdf(Omega_M | 4, diag_matrix(rep_vector(1,3)));
  target += inv_wishart_lpdf(Omega_F | 4, diag_matrix(rep_vector(1,3)));
  
  // Random-effects
  for(i in 1:n_M){ target += multi_normal_lpdf(bi_M[i,1:3] | rep_vector(0,3), Omega_M); }
  for(i in 1:n_F){ target += multi_normal_lpdf(bi_F[i,1:3] | rep_vector(0,3), Omega_F); }
  
  // Survival fixed effects
  target += normal_lpdf(beta_2 | 0, 10);
  target += normal_lpdf(gamma_2 | 0, 10);
 
  // Association parameters
  target += normal_lpdf(alphaB_2 | 0, 10);
  target += normal_lpdf(alphaG_2 | 0, 10);
  target += normal_lpdf(alphaD_2 | 0, 10);

  // Shape parameters (Weibull hazard)
  target += cauchy_lpdf(phi_2 | 0, 1);
  
}

"


# TWO-STAGE SPECIFICATION

BiExp <- "
// BI-EXPONENTIAL SPECIFICATION
functions {
  vector nonlinear_predictor(int[] ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int N;
  int n;
  vector[N] y;
  int<lower=1,upper=n> ID[N];
  vector[N] times;
  int<lower=1,upper=N> start[n];
  int<lower=1,upper=N> stop[n];
}
parameters {
  vector[3] theta;
  real<lower=0> sigma2;
  cov_matrix[3] Omega;
  matrix[n,3] bi;
}
model {
  // BI-EXPONENTIAL SPECIFICATION
  vector[N] nonlinpred = nonlinear_predictor(ID, times, theta, bi);
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y | nonlinpred, sqrt(sigma2));

  // LOG-PRIORS
  // Longitudinal fixed effects
  target += normal_lpdf(theta | 0, 10);
   
  // Residual error variance
  target += cauchy_lpdf(sigma2 | 0, 5);
 
  // Random-effects variance-covariance matrix
  target += inv_wishart_lpdf(Omega | 4, diag_matrix(rep_vector(1,3)));
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:3] | rep_vector(0,3), Omega); }
   
}

"


WeibPHCompRisk <- "
// BI-EXPONENTIAL SPECIFICATION
functions {
  vector nonlinear_predictor(int[] ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int n;
  int nbeta;
  matrix[n,nbeta] X;
  vector[n] Time;
  vector[n] status_1;
  vector[n] status_2;
  // Longitudinal data 
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  vector[N_M] y_M;
  vector[N_F] y_F;
  int<lower=1,upper=n_M> ID_M1[N_M];
  int<lower=1,upper=n> ID_M2[N_M];
  int<lower=1,upper=n_F> ID_F1[N_F];
  int<lower=1,upper=n> ID_F2[N_F];
  vector[N_M] times_M;
  vector[N_F] times_F;
  // Shared information 
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  matrix[3,3] Omega_M;
  matrix[3,3] Omega_F;
}
parameters {
  vector[nbeta] beta_1;
  vector[nbeta] beta_2;
  vector[2] alphaB_1;
  vector[2] alphaG_1;
  vector[2] alphaD_1;
  vector[2] alphaB_2;
  vector[2] alphaG_2;
  vector[2] alphaD_2;
  real gamma_1;
  real gamma_2;
  real<lower=0> phi_1;
  real<lower=0> phi_2;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
model {
  vector[n] loghaz_1;
  vector[n] loghaz_2;
  vector[n] cumHaz_1;
  vector[n] cumHaz_2;
  
  // BI-EXPONENTIAL SPECIFICATION
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y_M | nonlinpred_M, sqrt(sigma2_M));
  target += normal_lpdf(y_F | nonlinpred_F, sqrt(sigma2_F));
  
  // COMPETING RISKS LOG-LIKELIHOOD
  for(i in 1:n){
       // Hazard functions
       loghaz_1[i] = log(phi_1) + (phi_1-1)*log(Time[i]) + gamma_1 + X[i,]*beta_1 + 
                       lBi[i,]*alphaB_1 + lGi[i,]*alphaG_1 + lDi[i,]*alphaD_1;
       loghaz_2[i] = log(phi_2) + (phi_2-1)*log(Time[i]) + gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2;

       // Cumulative hazard functions
       cumHaz_1[i] = pow(Time[i], phi_1) * exp( gamma_1 + X[i,]*beta_1 + 
                       lBi[i,]*alphaB_1 + lGi[i,]*alphaG_1 + lDi[i,]*alphaD_1 );
       cumHaz_2[i] = pow(Time[i], phi_2) * exp( gamma_2 + X[i,]*beta_2 + 
                       lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2 );
       
       target += status_1[i]*loghaz_1[i] + status_2[i]*loghaz_2[i] - cumHaz_1[i] - cumHaz_2[i];
  }
  
  // LOG-PRIORS
  // Survival fixed effects
  target += normal_lpdf(beta_1 | 0, 10);  
  target += normal_lpdf(beta_2 | 0, 10);
  target += normal_lpdf(gamma_1 | 0, 10);  
  target += normal_lpdf(gamma_2 | 0, 10);
 
  // Association parameters
  target += normal_lpdf(alphaB_1 | 0, 10);
  target += normal_lpdf(alphaB_2 | 0, 10);
  target += normal_lpdf(alphaG_1 | 0, 10);
  target += normal_lpdf(alphaG_2 | 0, 10);
  target += normal_lpdf(alphaD_1 | 0, 10);
  target += normal_lpdf(alphaD_2 | 0, 10);

  // Shape parameters (Weibull hazard)
  target += cauchy_lpdf(phi_1 | 0, 1);
  target += cauchy_lpdf(phi_2 | 0, 1);
  
  // Individual fixed-effects
  for(i in 1:n_M){ target += multi_normal_lpdf(bi_M[i,1:3] | rep_vector(0,3), Omega_M); }
  for(i in 1:n_F){ target += multi_normal_lpdf(bi_F[i,1:3] | rep_vector(0,3), Omega_F); }
   
}

"


WeibPHSurv <- "
// BI-EXPONENTIAL SPECIFICATION
functions {
  vector nonlinear_predictor(int[] ID, vector times, vector theta, matrix bi){
     int N = num_elements(times);
     vector[N] logBi = theta[1] + bi[ID,1];      // log-baseline
     vector[N] Gi = exp(theta[2] + bi[ID,2]);    // growth rate
     vector[N] Di = exp(theta[3] + bi[ID,3]);    // decay rate

     vector[N] out = logBi + log( exp(rows_dot_product(Gi, times)) + exp(-rows_dot_product(Di, times)) - 1 );

     return out;
  }
}
data {
  int n;
  int nbeta;
  matrix[n,nbeta] X;
  vector[n] Time;
  vector[n] status_2;
  // Longitudinal data 
  int N_M;
  int N_F;
  int n_M;
  int n_F;
  vector[N_M] y_M;
  vector[N_F] y_F;
  int<lower=1,upper=n_M> ID_M1[N_M];
  int<lower=1,upper=n> ID_M2[N_M];
  int<lower=1,upper=n_F> ID_F1[N_F];
  int<lower=1,upper=n> ID_F2[N_F];
  vector[N_M] times_M;
  vector[N_F] times_F;
  // Shared information 
  vector[3] theta_M;
  vector[3] theta_F;
  real<lower=0> sigma2_M;
  real<lower=0> sigma2_F;
  matrix[3,3] Omega_M;
  matrix[3,3] Omega_F;
}
parameters {
  vector[nbeta] beta_2;
  vector[2] alphaB_2;
  vector[2] alphaG_2;
  vector[2] alphaD_2;
  real gamma_2;
  real<lower=0> phi_2;
  matrix[n_M,3] bi_M;
  matrix[n_F,3] bi_F;
}
transformed parameters {
  // log-baseline, log-growth and log-decay paramaters
  matrix[n,2] lBi;
  matrix[n,2] lGi;
  matrix[n,2] lDi;
  
  // M-spike
  lBi[,1] = rep_vector(theta_M[1],n);
  lBi[ID_M2,1] = theta_M[1] + bi_M[ID_M1,1];
  lGi[,1] = rep_vector(theta_M[2],n);
  lGi[ID_M2,1] = theta_M[2] + bi_M[ID_M1,2];
  lDi[,1] = rep_vector(theta_M[3],n);
  lDi[ID_M2,1] = theta_M[3] + bi_M[ID_M1,3];
  
  // FLC
  lBi[,2] = rep_vector(theta_F[1],n);
  lBi[ID_F2,2] = theta_F[1] + bi_F[ID_F1,1];
  lGi[,2] = rep_vector(theta_F[2],n);
  lGi[ID_F2,2] = theta_F[2] + bi_F[ID_F1,2];
  lDi[,2] = rep_vector(theta_F[3],n);
  lDi[ID_F2,2] = theta_F[3] + bi_F[ID_F1,3];
}
model {
  vector[n] loghaz_2;
  vector[n] cumHaz_2;
   
  // BI-EXPONENTIAL SPECIFICATION
  vector[N_M] nonlinpred_M = nonlinear_predictor(ID_M1, times_M, theta_M, bi_M);
  vector[N_F] nonlinpred_F = nonlinear_predictor(ID_F1, times_F, theta_F, bi_F);
   
  // LONGITUDINAL NORMAL LOG-LIKELIHOOD
  target += normal_lpdf(y_M | nonlinpred_M, sqrt(sigma2_M));
  target += normal_lpdf(y_F | nonlinpred_F, sqrt(sigma2_F));   
   
  // SURVIVAL LOG-LIKELIHOOD
  for(i in 1:n){
       // Hazard function
       loghaz_2[i] = log(phi_2) + (phi_2-1)*log(Time[i]) + gamma_2 + X[i,]*beta_2 + 
                    lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2;

       // Cumulative hazard functions
       cumHaz_2[i] = pow(Time[i], phi_2) * exp( gamma_2 + X[i,]*beta_2 + 
                    lBi[i,]*alphaB_2 + lGi[i,]*alphaG_2 + lDi[i,]*alphaD_2 );
       
       target += status_2[i]*loghaz_2[i] - cumHaz_2[i];
  }
  
  // LOG-PRIORS
  // Survival fixed effects
  target += normal_lpdf(beta_2 | 0, 10);
  target += normal_lpdf(gamma_2 | 0, 10);
 
  // Association parameters (current info)
  target += normal_lpdf(alphaB_2 | 0, 10);
  target += normal_lpdf(alphaG_2 | 0, 10);
  target += normal_lpdf(alphaD_2 | 0, 10);

  // Shape parameters (Weibull hazard)
  target += cauchy_lpdf(phi_2 | 0, 1);
  
  // Individual fixed-effects
  for(i in 1:n_M){ target += multi_normal_lpdf(bi_M[i,1:3] | rep_vector(0,3), Omega_M); }
  for(i in 1:n_F){ target += multi_normal_lpdf(bi_F[i,1:3] | rep_vector(0,3), Omega_F); }
   
}

"
