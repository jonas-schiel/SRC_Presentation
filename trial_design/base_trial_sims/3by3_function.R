#Function to simulate a 3+3 trial

##############################Trial Rules#######################################
#
#1. Doses administered to cohorts with a size of 3 patients
#2. Doses can never re-escalate after de-escalating
#3. Each dose can only be administered a total of 6 times (2 cohorts)
#4. Dose is disqualified when DLT greater than 1/3 or 1/6
#
################################################################################



##############################Parameters########################################
#
# doses = Dose levels that will be administered (entered as a vector, c())
#
# dlt_probabilities = Probability that each dose might result in adverse 
# reactions (vector that corresponds to doses vector)
#
################################################################################


ThreeThree_sim_func <- function(doses, dlt_probs) {
  n_j <- rep(0, length(doses))
  y_sum_j <- rep(0, length(doses))
  i <- 1
  desc <- 0
  co <- 0
  n <- 0
  
  ThreeThree_data <- data.frame(
    cohort = integer(),
    patient = integer(),
    dose = integer(),
    dlt = integer()
  )
  
  while (n_j[i] < 6 && desc == 0) {
    co <- co + 1
    trial <- rbinom(3, 1, dlt_probs[i])
    
    n_j[i] <- n_j[i] + length(trial)
    y_sum_j[i] <- y_sum_j[i] + sum(trial)
    
    row_1 <- data.frame(
      cohort = co,
      patient = n*3 + 1,
      dose = doses[i],
      dlt = trial[1]
    )
    
    row_2 <- data.frame(
      cohort = co,
      patient = n*3 + 2,
      dose = doses[i],
      dlt = trial[2]
    )
    
    row_3 <- data.frame(
      cohort = co,
      patient = n*3 + 3,
      dose = doses[i],
      dlt = trial[3]
    )
    
    ThreeThree_data <- rbind(ThreeThree_data, row_1)
    ThreeThree_data <- rbind(ThreeThree_data, row_2)
    ThreeThree_data <- rbind(ThreeThree_data, row_3)
    
    escalate <- ifelse(y_sum_j[i] == 0 && n_j[i] == 3 && i < length(doses) || (y_sum_j[i] / n_j[i]) <= 1/6 && n_j[i] == 6 && i < length(doses), 1, 0)
    stay <- ifelse(y_sum_j[i] == 1 && n_j[i] == 3 || 
                     (y_sum_j[i] / n_j[i]) <= 1/6 && n_j[i] == 6 && i == length(doses) ||
                     y_sum_j[i] == 0 && n_j[i] == 3 && i == length(doses), 1, 0)
    decrease <- ifelse(y_sum_j[i] > 1 && n_j[i] == 3 && i > 1 || y_sum_j[i] > 1 && n_j[i] == 6 && i > 1, 1, 0)
    terminate <- ifelse(y_sum_j[i] > 1 && n_j[i] == 3 && i == 1 || y_sum_j[i] > 1 && n_j[i] == 6 && i == 1, 1, 0)
    
    if (escalate == 1) {
      i <- i+1
    } else if (stay == 1) {
      i <- i
    } else if (decrease == 1) {
      i <- i - 1
    } else if (terminate == 1) {
      i <- 101
    } else {
      i <- 401
    }
    
    if (decrease == 1) {
      desc <- 1
    }
    
    n <- n + 1
    
    if (i == 101 || i == 401) {
      break
    }
  }
  
  MTD <- ifelse(terminate == 1, NA,
                ifelse(decrease == 1, doses[i],
                       doses[i+1]))
  
  if (i == 401) {
    cat("Coding Error")
  } else {
    results <- list(
      ThreeThree_data = ThreeThree_data,
      MTD = MTD
    )
    return(results)
  }
}
