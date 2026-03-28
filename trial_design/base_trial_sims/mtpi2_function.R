#Simulates an mTPI-2 trial (modified toxicity probability interval-2)

############################Rules for the trial#################################
#
#1. Trial cannot administer more than n_max patients 
#2. Trial will end once a 95% credible interval at dose_j (last administered 
#   dose) falls below a given width
#3. If the posterior probability of a DLT at dose_j exceeds 65%, the dose is
#   removed from the trial
#
################################################################################



#############################Function Parameters################################
# 
# p_t = Target toxicity
# 
# e1 = Error below p_t, determines how much error below the target toxicity is
# acceptable
# 
# e2 = Error above p_t, determines how much error above the target toxcity is
# acceptable
# 
# doses = Dose levels that will be administered (entered as a vector, c())
# 
# dlt_probabilities = Probability that each dose might result in adverse
# reactions (vector that corresponds to doses vector)
# 
# n_max = Maximum number of patients that can be administered to the trial
# 
# ci_end_width = Width of 95% Credible Interval at dose_j required to conclude
# the trial
#
################################################################################

mTPI2_sim_func <- function(p_t, e1, e2, doses, dlt_probs, n_max, ci_end_width) {
  buffer <- e1 + e2
  l <- p_t - e1
  u <- p_t + e2
  
  upper_vector <- seq(u, 1, by = buffer)
  if (upper_vector[length(upper_vector)] < 1) {
    upper_vector[length(upper_vector) + 1] <- 1
  }
  lower_vector <- seq(l, 0, by = -buffer)
  if (lower_vector[length(lower_vector)] > 0) {
    lower_vector[length(lower_vector) + 1] <- 0
  }
  
  n_j <- rep(0, length(doses))
  y_sum_j <- rep(0, length(doses))  
  post_means <- rep(NA, length(doses))
  i <- 1
  n <- 0
  ci_width <- 1
  
  mTPI2_trial_data <- data.frame(
    patient = integer(),
    dose = integer(),
    dlt = integer(),
    decision = character()
  )
  
  while (n < n_max && ci_width > ci_end_width) {
    n <- n + 1
    
    trial <- rbinom(1, 1, dlt_probs[i])
    n_j[i] <- n_j[i] + length(trial)
    y_sum_j[i] <- y_sum_j[i] + sum(trial)
    
    a_star <- 0.5 + y_sum_j[i]
    b_star <- 0.5 + n_j[i] - y_sum_j[i]
    post_means[i] <- a_star / (a_star + b_star)
    post_cred <- qbeta(c(0.025, 0.975), a_star, b_star)
    ci_width <- diff(post_cred)
    
    s_prob <- pbeta(u, a_star, b_star, lower.tail = TRUE) - pbeta(l, a_star, b_star, lower.tail = TRUE)
    
    k <- length(lower_vector)
    lower_post_vec <- rep(0, length(lower_vector)-1)
    lower_int_length_vec <- rep(0, length(lower_vector)-1)
    while (k > 1) {
      val <- pbeta(lower_vector[k-1], a_star, b_star, lower.tail = TRUE) - pbeta(lower_vector[k], a_star, b_star, lower.tail = TRUE)
      int_length <- lower_vector[k-1] - lower_vector[k]
      lower_post_vec[k] <- val
      lower_int_length_vec[k] <- int_length
      k <- k - 1
    }
    
    lower_int_length_vec <- lower_int_length_vec[2:length(lower_int_length_vec)]
    lower_post_vec <- lower_post_vec[2:length(lower_post_vec)]
    
    k <- 1
    upper_post_vec <- rep(0, length(upper_vector)-1)
    upper_int_length_vec <-  rep(0, length(upper_vector)-1)
    while (k < length(upper_vector)) {
      val <- pbeta(upper_vector[k + 1], a_star, b_star, lower.tail = TRUE) - pbeta(upper_vector[k], a_star, b_star, lower.tail = TRUE)
      int_length <- upper_vector[k + 1] - upper_vector[k]
      upper_post_vec[k] <- val
      upper_int_length_vec[k] <- int_length
      k <- k + 1
    }
    
    s_norm <- s_prob / buffer
    e_norm <- lower_post_vec / lower_int_length_vec
    d_norm <- upper_post_vec / upper_int_length_vec
    
    e_max <- max(e_norm)
    d_max <- max(d_norm)
    
    upm <- c(e_max, s_norm, d_max)
    result <- which.max(upm)
    result < as.integer(result)
    
    dec <- ifelse(result == 1 && i < length(doses), "escalate",
                  ifelse(result == 1 && i == length(doses) || result == 2 && i == length(doses), "stay (at max)",
                         ifelse(result == 2, "stay",
                                ifelse(result == 3 && i > 1, "de-escalate",
                                       "terminate trial"))))
    
    new_data <- data.frame(
      patient = n,
      DLT = sum(trial),
      dose = doses[i],
      decision = dec,
      posterior_mean = post_means[i]
    )
    
    mTPI2_trial_data <- rbind(mTPI2_trial_data, new_data)
    
    posterior_safety_check <- pbeta(p_t, a_star, b_star, lower.tail = FALSE)
    if (posterior_safety_check >= 0.65 && i > 1 && n_j[i] > 2) {
      doses <- doses[1:(i-1)]
    }
    if (posterior_safety_check >= 0.65 && i == 1 && n_j[i] > 2) {
      # cat("Trial Terminated. Dose 1 exceeds MTD.")
      stop <- 1
      break
    }
    
    i <- ifelse(result == 1 && i < length(doses), i + 1,
                ifelse(result == 1 && i == length(doses), i, 
                       ifelse(result == 2, i,
                              ifelse(result == 3 && i > 1, i - 1,
                                     ifelse(result == 3 && i == 1, 101,
                                            401)))))
    
    
    
    if (i == 101) {
      break
    }
    if (i == 401) {
      break
    }
  }
  
  if (i == 101) {
    # cat("Trial Suspended")
  }
  if (i == 401) {
    cat("Coding Error")
  }
  
  a_star_vec <- 0.5 + y_sum_j
  b_star_vec <- 0.5 + n_j - y_sum_j
  posterior_probability_DLT <- pbeta(p_t, a_star_vec, b_star_vec, lower.tail = FALSE)
  
  posterior_probability_DLT <- posterior_probability_DLT[1:length(doses)]
  
  posterior_df <- data.frame(
    dose = doses,
    posterior_probability = posterior_probability_DLT
  )
  
  if (i != 101) {
    MTD <- mTPI2_trial_data[nrow(mTPI2_trial_data), "dose"]
  } else {
    MTD <- NA
  }
  
  return(list(
    patient_data = mTPI2_trial_data,
    probability_data = posterior_df,
    MTD = MTD
  ))
}