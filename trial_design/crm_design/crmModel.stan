
data { //declaring what we will be putting into the model
//Trial data here:
//2. Total number of doses that can potntially be tested in trial
  int<lower=0> doses;
//3. Arrays to store number of observations and DLTs at each dose
  vector[doses] dose_level;
  array[doses] int<lower=0> n_at_dose;
  array[doses] int<lower=0> dlt_at_dose;


//Pseudo data for prior "Pseudo-Posterior"
//1. Total number of Pseudo observations
  int N_prior;
//2. Vector for each prior dose level
  vector[N_prior] dosing_prior;
//3. Vector for each prior dlt prob
  vector<lower=0,upper=1>[N_prior] dlt_prior;
//4. Array for effective sample size
  array[N_prior] real<lower=0> eff_n;
}

parameters {//outlining prior parameters
  real<upper=0> beta0;
  real<lower=0> beta1;
}

model {
  //give distrbutions to prior parameters
  beta0 ~ normal(-5, 2.5);
  beta1 ~ gamma(2, 100);
  
  //Weight the Pseudo data to the model
  for (i in 1:N_prior) {
    //convert log-odds to probabilty at each prior dose
    real prob = inv_logit(beta0 + beta1 * dosing_prior[i]);
    //multiply effective sample size by log-likelihood
    target += eff_n[i] * (dlt_prior[i]*log(prob) + ((1-dlt_prior[i])*log1m(prob)));
  }
  //handle the trial data
  dlt_at_dose ~ binomial_logit(n_at_dose, beta0 + beta1 * dose_level);
}

