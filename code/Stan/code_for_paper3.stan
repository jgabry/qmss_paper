data {
  // dimensions 
  int<lower=1> N ; # number of observations 
  int<lower=1> C ; # number of congresses (time periods)
  
  // arrays of variables 
  int<lower=1,upper=C> cong[N] ;  # map between observations & congresses
  int<lower=1,upper=56> nVotes[N] ;  # total votes
  int<lower=0,upper=55> nWins[N] ;  # majority party victories
  real lvRatio[N] ;  # log(ratio of avg vote-shares)
  
  // inverse of penalty matrix 
  matrix[C,C] Pinverse ;  # for GMRF prior
}
transformed data {
  real<lower=0> phi_scale ;  # scale for prior on phi
  real<lower=0> phi_loc ;  # location for prior on phi
  real<lower=0> tau_scale ;  # scale for priors on taus
  real<lower=0> tau_loc ;  # location for priors on taus
  matrix[C,C]   cholPinverse ;  # Cholesky decomposition 
  
  phi_loc <- 0 ;
  phi_scale <- 10 ;
  tau_loc <- 0 ;
  tau_scale <- 1 ;
  cholPinverse <- cholesky_decompose(Pinverse) ;
}
parameters {
  vector[C] lambda_noise ;      
  vector[C] rho_noise ;
  real<lower=0> phi_noise ;    
  real<lower=0> tau_sq_lambda ;  
  real<lower=0> tau_sq_rho ;
}
transformed parameters {
  real<lower=0> phi ;
  vector[C] lambda ;      
  vector[C] rho ;
  
  // inverse CDF method, standard normal --> half-Cauchy
  phi <- phi_loc + phi_scale * tan(pi() * (Phi_approx(phi_noise) - 0.5)) ;
  
  // non-centered parameterization
  lambda  <- (tau_sq_lambda * cholPinverse) * lambda_noise ;
  rho  <- (tau_sq_rho * cholPinverse) * rho_noise ;
}
model {
  // local/temporary variables
  real logLik ;  # log likelihood
  real logPrior ;  # log prior
  vector<lower=0>[N] alpha ;  # shape1 for Beta-Binomial 
  vector<lower=0>[N] beta ;  # shape2 for Beta-Binomial

  // log priors
  logPrior <- (
    normal_log(phi_noise, 0, 1) +
    normal_log(tau_sq_lambda, tau_loc, tau_scale) +
    normal_log(tau_sq_rho, tau_loc, tau_scale) +
    normal_log(lambda_noise, 0, 1) + 
    normal_log(rho_noise, 0, 1) 
  ) ;
  
  // log likelihood
  for (n in 1:N) {
    real eta_n ;
    real theta_n ;
    eta_n <- lambda[cong[n]] + rho[cong[n]] * lvRatio[n] ;
    theta_n   <- inv_logit(eta_n) ;    
    alpha[n] <- theta_n * phi ;
    beta[n]  <- (1 - theta_n) * phi ;
  } 
  logLik <- beta_binomial_log(nWins, nVotes, alpha, beta) ; # vectorized 
  
  // log posterior (up to proportion)
  increment_log_prob(logPrior + logLik) ; 
}
generated quantities {
  real bias[C] ;  # estimated bias toward majority party
  int y_rep[N] ;  # simulated data from posterior predictive dist.
  real theta[N] ;  
  
  for (c in 1:C) 
    bias[c] <- inv_logit(lambda[c]) - 0.5 ;
  
  for (n in 1:N) {
    real eta_n ;
    real alpha_n ;
    real beta_n ;
    eta_n <- lambda[cong[n]] + rho[cong[n]] * lvRatio[n] ;
    theta[n] <- inv_logit(eta_n) ;    
    alpha_n <- phi * theta[n] ;
    beta_n  <- phi * (1 - theta[n]) ;
    y_rep[n] <- beta_binomial_rng(nVotes[n], alpha_n, beta_n) ;
  }
}
