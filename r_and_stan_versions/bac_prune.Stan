// LD pruning homemade

data {
  int<lower=1> N;
  int<lower=1> V;
  matrix[N,V] gt;
}

parameters {
  real alpha;
  vector[V] beta_SNPs;
  real<lower=0> sigma;
}

model {
  vector[N] mu = beta_SNPs * gt + alpha;
  
  alpha ~ normal(0, 5);
  beta_SNPs ~ normal(0, 5);
  
  target += normal_lpdf(gt | mu, sigma);
  
}

