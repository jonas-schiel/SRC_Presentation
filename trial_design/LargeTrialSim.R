# Function to simulate trials at scale and collect data

#############################Function Parameters################################
# 
# doses = Dose levels that will be administered (entered as a vector, c())
# 
# dlt_probabilities = Probability that each dose might result in adverse
# reactions (vector that corresponds to doses vector)
#
# p_t = Target toxicity
# 
# e1 = Error below p_t, determines how much error below the target toxicity is
# acceptable
# 
# e2 = Error above p_t, determines how much error above the target toxcity is
# acceptable
# 
# ci_end_width = Width of 95% Credible Interval at dose_j required to conclude
# the trial
#
# n_max = Maximum number of patients that can be administered to the trial
#
# sims = total number of simulations per each model
#
################################################################################

get_data <- function(doses, dlt_probs, p_t, e1, e2, ci_end_width, n_max, sims) {
  
  source("base_trial_sims/3by3_function.R")
  source("base_trial_sims/mtpi_function.R")
  source("base_trial_sims/mtpi2_function.R")
  
  ThreeThree_data <- data.frame(
    patients = integer(),
    trial_dlts = integer(),
    trial_toxicity = numeric(),
    MTD = integer()
  )
  
  mtpi_data <- data.frame(
    patients = integer(),
    trial_dlts = integer(),
    trial_toxicity = numeric(),
    MTD = integer()
  )
  
  mtpi2_data <- data.frame(
    patients = integer(),
    trial_dlts = integer(),
    trial_toxicity = numeric(),
    MTD = integer()
  )
  
  for (i in 1:sims) {
    threethree_trial <- ThreeThree_sim_func(doses, dlt_probs)
    
    mtpi_trial <- mTPI_sim_func(p_t, err = e1, doses = doses, dlt_probs = dlt_probs, n_max = n_max, ci_end_width = ci_end_width)
    
    mtpi2_trial <- mTPI2_sim_func(p_t = p_t, e1 = e1, e2 = e2, doses = doses, dlt_probs = dlt_probs, n_max = n_max, ci_end_width = ci_end_width)
    
    dat_1 <- threethree_trial$ThreeThree_data
    threethree_pat <- nrow(dat_1)
    threethree_dlt <- sum(dat_1$dlt)
    threethree_MTD <- threethree_trial$MTD
    threethree_tt <- threethree_dlt / threethree_pat
    
    new_three_data <- data.frame(
      patients = threethree_pat,
      trial_dlts = threethree_dlt,
      trial_toxicity = threethree_tt,
      MTD = threethree_MTD
    )
    
    dat_2 <- mtpi_trial$patient_data
    mtpi_pat <- nrow(dat_2)
    mtpi_dlt <- sum(dat_2$DLT)
    mtpi_MTD <- mtpi_trial$MTD
    mtpi_tt <- mtpi_dlt / mtpi_pat
    
    new_mtpi_data <- data.frame(
      patients = mtpi_pat,
      trial_dlts = mtpi_dlt,
      trial_toxicity = mtpi_tt,
      MTD = mtpi_MTD
    )
    
    dat_3 <- mtpi2_trial$patient_data
    mtpi2_pat <- nrow(dat_3)
    mtpi2_dlt <- sum(dat_3$DLT)
    mtpi2_MTD <- mtpi2_trial$MTD
    mtpi2_tt <- mtpi2_dlt / mtpi2_pat
    
    new_mtpi2_data <- data.frame(
      patients = mtpi2_pat,
      trial_dlts = mtpi2_dlt,
      trial_toxicity = mtpi2_tt,
      MTD = mtpi2_MTD
    )
    
    
    ThreeThree_data <- rbind(ThreeThree_data, new_three_data)
    mtpi_data <- rbind(mtpi_data, new_mtpi_data)
    mtpi2_data <- rbind(mtpi2_data, new_mtpi2_data)
  }
  return(list(
    threethree = ThreeThree_data,
    mtpi = mtpi_data,
    mtpi2 = mtpi2_data
  ))
} 