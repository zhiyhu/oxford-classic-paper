/*A simple example of an hierarchical model*/
// https://biologyforfun.wordpress.com/2016/11/10/hierarchical-models-with-rstan-part-1/
data {
  int<lower=1> N; //the number of observations
  int<lower=1> J; //the number of groups
  int<lower=1,upper=J> id[N]; //vector of group indeces
  vector[N] X; //the model matrix
  vector[N] y; //the response variable
}
parameters {
  real gamma; //population-level regression coefficients
  real<lower=0> tau; //the standard deviation of the regression coefficients

  real beta[J]; //matrix of group-level regression coefficients
  real<lower=0> sigma; //standard deviation of the individual observations
}
model {
  vector[N] mu; //linear predictor
  //priors
  gamma ~ normal(0,5); //weakly informative priors on the regression coefficients
  tau ~ cauchy(0,2.5); //weakly informative priors, see section 6.9 in STAN user guide
  sigma ~ gamma(2,0.1); //weakly informative priors, see section 6.9 in STAN user guide
  
  for(j in 1:J){
   beta[j] ~ normal(gamma,tau); //fill the matrix of group-level regression coefficients 
  }
  
  for(n in 1:N){
    mu[n] = X[n] * beta[id[n]]; //compute the linear predictor using relevant group-level regression coefficients 
  }

  //likelihood
  y ~ normal(mu,sigma);
}
