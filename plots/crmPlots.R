# Code used to plot the CRM model plots and trial results

source("../trial_design/crm_design/crmSimulation.R")


##################Get Prior Distribution DataFrame##############################

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

stan_prior <- extract(prior.model)

prior_func <- sapply(1:360, function(s) {
  plogis(stan_prior$beta0 + stan_prior$beta1*s)
})

stan_mean <- apply(prior_func, 2, mean)
stan_upper <- apply(prior_func, 2, quantile, 0.95)
stan_lower <- apply(prior_func, 2, quantile, 0.05)

stan_df <- data.frame(
  mean = stan_mean,
  upper = stan_upper,
  lower = stan_lower
)

################################################################################


####################Simulate Trial Under Various Toxicities#####################

options(mc.cores = parallel::detectCores())
high <- crm_data_sim(stan_file = here("trial_design/crm_design/crmModel.stan"), burn_ins = 2000, iters = 5000, chains = 3,
                     t_beta0 = -3, t_beta1 = 0.07, dose_level = c(10, 20, 40, 80, 160, 320),
                     end_ci_width = 0.30, ci_confidence = 0.90, p_t = 0.3,
                     n_prior = 6, dlt_priors = c(0.05, 0.075, 0.10, 0.18, 0.30, 0.50),
                     eff_samp_priors = c(1, 0.75, 0.5, 0.25, 0.125, 0.0625))
moderate_high <- crm_data_sim(stan_file = here("trial_design/crm_design/crmModel.stan"), burn_ins = 2000, iters = 5000, chains = 3,
                              t_beta0 = -4.5, t_beta1 = 0.035, dose_level = c(10, 20, 40, 80, 160, 320),
                              end_ci_width = 0.30, ci_confidence = 0.90, p_t = 0.3,
                              n_prior = 6, dlt_priors = c(0.05, 0.075, 0.10, 0.18, 0.30, 0.50),
                              eff_samp_priors = c(1, 0.75, 0.5, 0.25, 0.125, 0.0625))

moderate_low <- crm_data_sim(stan_file = here("trial_design/crm_design/crmModel.stan"), burn_ins = 2000, iters = 5000, chains = 3,
                             t_beta0 = -5, t_beta1 = 0.025, dose_level = c(10, 20, 40, 80, 160, 320),
                             end_ci_width = 0.3, ci_confidence = 0.90, p_t = 0.3,
                             n_prior = 6, dlt_priors = c(0.05, 0.075, 0.10, 0.18, 0.30, 0.50),
                             eff_samp_priors = c(1, 0.75, 0.5, 0.25, 0.125, 0.0625))
low <- crm_data_sim(stan_file = here("trial_design/crm_design/crmModel.stan"), burn_ins = 2000, iters = 5000, chains = 3,
                    t_beta0 = -5, t_beta1 = 0.012, dose_level = c(10, 20, 40, 80, 160, 320),
                    end_ci_width = 0.3, ci_confidence = 0.90, p_t = 0.3,
                    n_prior = 6, dlt_priors =c(0.05, 0.075, 0.10, 0.18, 0.30, 0.50),
                    eff_samp_priors = c(1, 0.75, 0.5, 0.25, 0.125, 0.0625))

################################################################################



##################Extract Posteriors and Get Data For Curves####################

high_post <- high$posterior_data
moderate_high_post <- moderate_high$posterior_data
moderate_low_post <- moderate_low$posterior_data
low_post <- low$posterior_data


high_func <- sapply(1:360, function(s) {
  plogis(high_post$beta0 + high_post$beta1 * s)
})
moderate_high_func <- sapply(1:360, function(s) {
  plogis(moderate_high_post$beta0 + moderate_high_post$beta1 * s)
})
moderate_low_func <- sapply(1:360, function(s) {
  plogis(moderate_low_post$beta0 + moderate_low_post$beta1 * s)
})
low_func <- sapply(1:360, function(s) {
  plogis(low_post$beta0 + low_post$beta1 * s)
})

high_post_mean <- apply(high_func, 2, mean)
modhigh_post_mean <- apply(moderate_high_func, 2, mean)
modlow_post_mean <- apply(moderate_low_func, 2, mean)
low_post_mean <- apply(low_func, 2, mean)

high_upper <- apply(high_func, 2, quantile, 0.95)
modhigh_upper <- apply(moderate_high_func, 2, quantile, 0.95)
modlow_upper <- apply(moderate_low_func, 2, quantile, 0.95)
low_upper <- apply(low_func, 2, quantile, 0.95)

high_lower <- apply(high_func, 2, quantile, 0.05)
modhigh_lower <- apply(moderate_high_func, 2, quantile, 0.05)
modlow_lower <- apply(moderate_low_func, 2, quantile, 0.05)
low_lower <- apply(low_func, 2, quantile, 0.05)

high_tox_df <- data.frame(
  dose = 1:360,
  post_mean = high_post_mean,
  post_upper = high_upper,
  post_lower = high_lower,
  prior_mean = prior_mean,
  prior_upper = prior_upper,
  prior_lower = prior_lower
)



library(ggplot2)

################################################################################


###################Plotting under High Toxicity Assumption######################

high_tox_trial <- high$CRM_data

curve_plot_high <- ggplot(data = high_tox_df, aes(x = dose))+
  geom_line(aes(y = post_mean),col = "red")+
  scale_x_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(aes(y = post_upper),col = "red", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = post_lower),col = "red", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = prior_mean), col = "steelblue")+
  geom_line(aes(y = prior_upper), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_line(aes(y = prior_lower), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_hline(yintercept = 0.35, color = "grey", linetype = "dotted")+
  labs(title = "Informative Prior vs. Trial Posterior", subtitle = "When Treatment Toxicity is High", x = "Dose (mg)", y = "DLT Probability")+
  theme_classic()


trial_plot_high <- ggplot(data = high_tox_trial, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(dlt)))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30))+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  labs(title = "Simulated CRM Trial", subtitle = "When Treatment Toxicity is High", y = "Dose (mg)", x = "Patient")+
  theme_classic()


modhigh_tox_df <- data.frame(
  dose = 1:360,
  post_mean = modhigh_post_mean,
  post_upper = modhigh_upper,
  post_lower = modhigh_lower,
  prior_mean = prior_mean,
  prior_upper = prior_upper,
  prior_lower = prior_lower
)

################################################################################


###############Plotting under Moderate High Toxicity Assumption#################

modhigh_tox_trial <- moderate_high$CRM_data

curve_plot_modhigh <- ggplot(data = modhigh_tox_df, aes(x = dose))+
  geom_line(aes(y = post_mean),col = "orange")+
  scale_x_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(aes(y = post_upper),col = "orange", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = post_lower),col = "orange", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = prior_mean), col = "steelblue")+
  geom_line(aes(y = prior_upper), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_line(aes(y = prior_lower), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_hline(yintercept = 0.35, color = "grey", linetype = "dotted")+
  labs(title = "Informative Prior vs. Trial Posterior", subtitle = "When Treatment Toxicity is Moderate-High", x = "Dose (mg)", y = "DLT Probability")+
  theme_classic()


trial_plot_modhigh <- ggplot(data = modhigh_tox_trial, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(dlt)))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30))+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  labs(title = "Simulated CRM Trial", subtitle = "When Treatment Toxicity is Moderate-High", y = "Dose (mg)", x = "Patient")+
  theme_classic()

################################################################################


################Plotting under Moderate Low Toxicity Assumption#################

modlow_tox_df <- data.frame(
  dose = 1:360,
  post_mean = modlow_post_mean,
  post_upper = modlow_upper,
  post_lower = modlow_lower,
  prior_mean = prior_mean,
  prior_upper = prior_upper,
  prior_lower = prior_lower
)

modlow_tox_trial <- moderate_low$CRM_data

library(ggplot2)

curve_plot_modlow <- ggplot(data = modlow_tox_df, aes(x = dose))+
  geom_line(aes(y = post_mean),col = "blue")+
  scale_x_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(aes(y = post_upper),col = "blue", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = post_lower),col = "blue", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = prior_mean), col = "steelblue")+
  geom_line(aes(y = prior_upper), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_line(aes(y = prior_lower), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_hline(yintercept = 0.35, color = "grey", linetype = "dotted")+
  labs(title = "Informative Prior vs. Trial Posterior", subtitle = "When Treatment Toxicity is Moderate-Low", x = "Dose (mg)", y = "DLT Probability")+
  theme_classic()


trial_plot_modlow <- ggplot(data = modlow_tox_trial, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(dlt)))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30))+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  labs(title = "Simulated CRM Trial", subtitle = "When Treatment Toxicity is Moderate-Low", y = "Dose (mg)", x = "Patient")+
  theme_classic()

################################################################################


####################Plotting under Low Toxicity Assumption######################

low_tox_df <- data.frame(
  dose = 1:360,
  post_mean = low_post_mean,
  post_upper = low_upper,
  post_lower = low_lower,
  prior_mean = prior_mean,
  prior_upper = prior_upper,
  prior_lower = prior_lower
)

low_tox_trial <- low$CRM_data

library(ggplot2)

curve_plot_low <- ggplot(data = low_tox_df, aes(x = dose))+
  geom_line(aes(y = post_mean),col = "purple")+
  scale_x_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(aes(y = post_upper),col = "purple", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = post_lower),col = "purple", linetype = "dashed", alpha = 0.4)+
  geom_line(aes(y = prior_mean), col = "steelblue")+
  geom_line(aes(y = prior_upper), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_line(aes(y = prior_lower), col = "steelblue", linetype = "dashed", alpha = 0.25)+
  geom_hline(yintercept = 0.35, color = "grey", linetype = "dotted")+
  labs(title = "Informative Prior vs. Trial Posterior", subtitle = "When Treatment Toxicity is Low", x = "Dose (mg)", y = "DLT Probability")+
  theme_classic()


trial_plot_low <- ggplot(data = low_tox_trial, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(dlt)))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30))+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  labs(title = "Simulated CRM Trial", subtitle = "When Treatment Toxicity is Low", y = "Dose (mg)", x = "Patient")+
  theme_classic()

################################################################################

curve_plot_high
curve_plot_modhigh
curve_plot_modlow
curve_plot_low

trial_plot_high
trial_plot_modhigh
trial_plot_modlow
trial_plot_low