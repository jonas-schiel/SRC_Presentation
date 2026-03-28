#Simulates an mTPI trial (modified toxicity probability interval)

############################Rules for the trial#################################
#
#1. Trial cannot administer more than n_max patients 
#2. Trial will end once an 95% credible interval at dose_j (last administered 
#   dose) falls below a given width
#3. If the posterior probability of a DLT at dose_j exceeds 65%, the dose is
#   removed from the trial
#
################################################################################



#############################Function Parameters################################
#
# p_t = Target toxicity
#
# err = Error, determines how large the stay interval will be. Ex. err=0.05 will 
# result in an interval width of 0.10
#
# doses = Dose levels that will be administered (entered as a vector, c())
#
# dlt_probabilities = Probability that each dose might result in adverse 
# reactions (vector that corresponds to doses vector)
#
# n_max = Maximum number of patients that can be administered to the trial
#
# ci_end_width = Width of 90% Credible Interval at dose_j required to conclude 
# the trial
#
################################################################################

mTPI_sim_func <- function(p_t, err, doses, dlt_probs, n_max, ci_end_width) {
  l <- p_t - err
  u <- p_t + err
  
  mTPI_trial_data <- data.frame(
    patient = integer(),
    DLT = integer(),
    dose = integer(),
    decision = character(),
    posterior_mean = numeric()
  )
  
  i <- 1
  n <- 0
  n_j <- rep(0, length(doses))
  y_sum_j <- rep(0, length(doses))
  post_means <- rep(NA, length(doses))
  stop <- 0
  ci_width <- 1
  
  while (n < n_max && ci_width > ci_end_width) {
    n <- n + 1
    
    trial <- rbinom(1, 1, dlt_probs[i])
    
    n_j[i] <- n_j[i] + length(trial)
    y_sum_j[i] <- y_sum_j[i] + sum(trial)
    
    a_star <- 0.5 + y_sum_j[i]
    b_star <- 0.5 + n_j[i] - y_sum_j[i]
    post_means[i] <- a_star / (a_star + b_star)
    
    p_e <- pbeta(l, a_star, b_star, lower.tail = TRUE)
    p_s <- pbeta(u, a_star, b_star, lower.tail = TRUE) - pbeta(l, a_star, b_star, lower.tail = TRUE)
    p_d <- pbeta(u, a_star, b_star, lower.tail = FALSE)
    
    e <- p_e / l
    s <- p_s / (u-l)
    d <- p_d / (1-u)
    
    upm <- c(e, s, d)
    result <- which.max(upm)
    result <- as.integer(result)
    
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
    
    mTPI_trial_data <- rbind(mTPI_trial_data, new_data)
    
    post_cred <- qbeta(c(0.025, 0.975), a_star, b_star)
    ci_width <- diff(post_cred)
    
    posterior_safety_check <- pbeta(p_t, a_star, b_star, lower.tail = FALSE)
    if (posterior_safety_check >= 0.65 && i > 1 && n_j[i] > 2) {
      doses <- doses[1:(i-1)]
      n_j <- n_j[1:(i-1)]
      y_sum_j <- y_sum_j[1:(i-1)]
      post_means <- post_means[1:(i-1)]
      i <- i - 1 
      next
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
  
  a_star_vec <- 0.5 + y_sum_j[i]
  b_star_vec <- 0.5 + n_j[i] - y_sum_j[i]
  posterior_probability_no_DLT <- pbeta(p_t, a_star_vec, b_star_vec, lower.tail = TRUE)
  
  posterior_probability_no_DLT <- posterior_probability_no_DLT[1:length(doses)]
  
  posterior_df <- data.frame(
    dose = doses,
    posterior_probability = posterior_probability_no_DLT
  )
  
  if (i != 101) {
    MTD <- mTPI_trial_data[nrow(mTPI_trial_data), "dose"]
  } else {
    MTD <- NA
  }
  
  return(list(
    patient_data = mTPI_trial_data,
    probability_data = posterior_df,
    MTD = MTD
  ))
  
}