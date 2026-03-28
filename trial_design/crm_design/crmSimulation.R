# R code used to simulate CRM trials using stan sampler

library(here)

#######Develop and Evaluate Skeletal Prior For the Dose Response Curve##########

x_vals <- seq(1, 360, by = 1)

beta0 <- rnorm(5000, mean = -5, 2.5)
beta1 <- rgamma(5000, shape = 2, rate = 100)

prior <- sapply(x_vals, function(s) {
  plogis(beta0 + beta1*s)
})

prior_mean <- apply(prior, 2, mean)
prior_upper <- apply(prior, 2, quantile, 0.95)
prior_lower <- apply(prior, 2, quantile, 0.05)

plot(x_vals, prior_mean, col = "black", type = "l", ylim = c(0, 1))
lines(prior_upper, col = "red", lty = 2)
lines(prior_lower, col = "red", lty = 2)
abline(h = 0.3, lty = 3, col = "grey")

################################################################################


#######Run Stan model with Prior Elicited Above and Evalute Diagnostics##########

library(rstan)

prior_data <- list(
  N = 0,
  doses = 0,
  dose_level = numeric(0),
  n_at_dose = integer(0),
  dlt_at_dose = integer(0),
  
  N_prior = 6,
  dosing_prior = c(10, 20, 40, 80, 160, 320),
  dlt_prior = c(0.05, 0.075, 0.10, 0.18, 0.30, 0.50),
  eff_n = c(1, 0.75, 0.5, 0.25, 0.125, 0.0625)
)

prior.model <- stan(
  file = here("trial_design/crm_design/crmModel.stan"),
  data = prior_data,
  iter = 5000,
  warmup = 2000,
  chains = 3,
  init = function() list(beta0 = -5, beta1 = 0.02),
  refresh = 0
)

prior.model
traceplot(prior.model, pars = "beta0")
traceplot(prior.model, pars = "beta1")

################################################################################


####Determine the Different Trial Toxicities Under Which to Test the Model######

x_val <- seq(1, 360, by = 1)
beta0_samp1 <- -3
beta1_samp1 <- 0.07
beta0_samp2 <- -4.5
beta1_samp2 <- 0.035
beta0_samp3 <- -5
beta1_samp3 <- 0.025
beta0_samp4 <- -5
beta1_samp4 <- 0.012

dist1 <- sapply(x_val, function(i) {
  plogis(beta0_samp1 + beta1_samp1 * i)
})
dist2 <- sapply(x_val, function(i) {
  plogis(beta0_samp2 + beta1_samp2 * i)
})
dist3 <- sapply(x_val, function(i) {
  plogis(beta0_samp3 + beta1_samp3 * i)
})
dist4 <- sapply(x_val, function(i) {
  plogis(beta0_samp4 + beta1_samp4 * i)
})

plot(x_val, dist1, type = "l", col = "red")
lines(dist2, col = "orange")
lines(dist3, col = "blue")
lines(dist4, col = "purple")
abline(h = 0.3, lty = 3, col = "steelblue")

################################################################################


#################Function Used to Simulate CRM Trials###########################

#############################Trial Rules########################################
#
# 1. Dose must be tested at least twice before the next dose can be considered
#    eligible
# 2. The expected value of the posterior distribution of the model at dose j
#    must be below the trial target toxicity
# 3. The proportion of the posterior distribution above the target toxicity must
#    be less than 0.35 in order to escalate dose
# 4. If the proportion of the posterior distribution above the target toxicity 
#    exceeds 0.60 the dose must be de-escalated
#
################################################################################


#############################Function Parameters################################
# 
# stan_file = .stan file to be used in the trial
#
# burn_ins = warmup used in sampling
#
# iters = iterations used in sampling
#
# chains = chains used in sampling
#
# t_beta0 = assumed true value of beta0 used to generate dlt probabilities
#
# t_beta1 = assumed true value of beta1 used to generate dlt_probabilites
#
# dose_level = doses used in the trial (in the form of a vector, c())
#
# end_ci_width = Credible Interval width required to end the trial
#
# ci_confidence = Credible Interval certainty
#
# p_t = target toxicity
#
# n_prior = quantity of prior doses used in pseudo data (should match length of
#            dlt_priors and eff_samp_priors)
#
# dlt_priors = expectations of treatment toxicity at each dose (vector, c(),
#           corresponding with dose_level)
#
# eff_samp_priors = effective sample sizes used in pseudo prior (vector, c(),
#           corresponding with dose_level)
#
################################################################################
crm_data_sim <- function(stan_file, burn_ins, iters, chains, t_beta0, t_beta1,
                         dose_level, end_ci_width, ci_confidence,
                         p_t, n_prior, dlt_priors, eff_samp_priors) {
  library(rstan)
  
  #1. Get stan model with Pseudo data
  pseudo_data <- list(
    N = 0,
    doses = 0,
    dose_level = numeric(0),
    n_at_dose = integer(0),
    dlt_at_dose = integer(0),
    #Pseudo Data
    N_prior = n_prior,
    dosing_prior = dose_level,
    dlt_prior = dlt_priors,
    eff_n = eff_samp_priors
  )
  
  options(mc.cores = parallel::detectCores())
  compiled_model <- stan_model(file = stan_file)
  stan.model <- sampling(
    compiled_model,
    data = pseudo_data,
    warmup = burn_ins,
    iter = iters,
    chains = chains,
    init = function() list(beta0 = -6, beta1 = 0.04545),
    refresh = 0
  )
  
  #2. Initialize Values
  i <- 1
  N <- 0
  dlts_j <- rep(0, length(dose_level))
  n_j <- rep(0, length(dose_level))
  ci_width <- 1
  doses <- length(dose_level)
  valid_doses <- length(dose_level)
  
  #3. Empty DF
  CRM_data <- data.frame(
    patient = integer(),
    dose = integer(),
    dlt = integer(),
    CI_width = numeric(),
    DLT_prob = numeric()
  )
  
  #4. Sim Data
  while ((ci_width > end_ci_width && N < 40) || n_j[i] < 2) {
    #Update N, get "true" DLT prob at dose, sim one trial
    N <- N+1
    dlt_prob <- plogis(t_beta0 + t_beta1*dose_level[i])
    trial <- rbinom(1, 1, prob = dlt_prob)
    
    #update dlt and n vectors
    n_j[i] <- n_j[i] + length(trial)
    dlts_j[i] <- dlts_j[i] + sum(trial)
    
    #5. Update STAN model
    new_data <- list(
      N = N,
      doses = doses,
      dose_level = dose_level,
      n_at_dose = n_j,
      dlt_at_dose = dlts_j,
      #Pseudo Data
      N_prior = n_prior,
      dosing_prior = dose_level,
      dlt_prior = dlt_priors,
      eff_n = eff_samp_priors
    )
    
    options(mc.cores = parallel::detectCores())
    new.stan <- sampling(
      compiled_model,
      data = new_data,
      warmup = burn_ins,
      iter = iters,
      chains = chains,
      init = function() list(beta0 = -6, beta1 = 0.04545),
      refresh = 0
    )
    
    #6. Extract the posterior predictive at dose
    posterior <- extract(new.stan)
    
    #posterior predictive mean at dose j
    post_pred_j <- plogis(posterior$beta0 + posterior$beta1 * dose_level[i])
    
    #posterior function
    max_dose <- which.max(dose_level)
    post_func <- sapply(1:dose_level[max_dose], function(s) {
      plogis(posterior$beta0 + posterior$beta1 * s)
    })
    
    #find CI width at target toxicity or last dose
    buffer <- (1-ci_confidence)/2
    post_mean <- apply(post_func, 2, mean)
    post_upper <- apply(post_func, 2, quantile, (1-buffer))
    post_lower <- apply(post_func, 2, quantile, buffer)
    
    # target_index <- which.min(abs(post_mean-p_t))
    # ci_width <- post_upper[target_index] - post_lower[target_index]
    current_dose <- dose_level[i]
    ci_width <- post_upper[current_dose] - post_lower[current_dose]
    
    #7. Update trial data
    new_data <- data.frame(
      patient = N,
      dose = dose_level[i],
      dlt = sum(trial),
      CI_width = ci_width,
      DLT_prob = mean(post_pred_j)
    )
    
    CRM_data <- rbind(CRM_data, new_data)
    
    #8. Assess rules for escalation
    #rule 1: at least 2 samples at dose
    if (n_j[i] < 2) {
      rule_1 <- 0
    }
    else {
      rule_1 <- 1
    }
    #rule 2: E(DLT_j) < target_toxicity]
    dlt_j <- mean(post_pred_j)
    if (dlt_j > p_t) {
      rule_2 <- 0
    }
    else {
      rule_2 <- 1
    }
    #rule 3: Pr(DLT_j > target toxicity) < 35%
    dlt_j_above_tox <- mean(post_pred_j > p_t)
    if (dlt_j_above_tox >= 0.35) {
      rule_3 <- 0
    }
    else {
      rule_3 <- 1
    }
    #rule 4: Pr(DLT_j > target toxicity) < 60% (de-escalation rule)
    if (dlt_j_above_tox < 0.60) {
      rule_4 <- 0
    }
    else {
      rule_4 <- 1
    }
    
    #9. Update i according to rules
    i <- ifelse(i < valid_doses && rule_1 == 1 && rule_2 == 1 && rule_3 == 1, i+1,
                i)
    i <- ifelse(i > 1 && rule_4 == 1, i - 1,
                ifelse(i == 1 && rule_4 == 1, 101,
                       i))
    
    valid_doses <- ifelse(rule_4 == 1, valid_doses - (valid_doses - i),
                          valid_doses)
    
    if (i == 101) {
      break
    }
    
  }
  
  #10 return dataframe 
  if (i == 101) {
    CRM_data[nrow(CRM_data), "dose"] <- NA
    prior_data <- extract(stan.model)
    post_data <- extract(new.stan)
    return(list(
      CRM_data = CRM_data,
      prior_data = prior_data,
      posterior_data = post_data,
      model = stan.model))
  } else {
    prior_data <- extract(stan.model)
    post_data <- extract(new.stan)
    return(list(
      CRM_data = CRM_data,
      prior_data = prior_data,
      posterior_data = post_data,
      model = stan.model))
  }
}

################################################################################



