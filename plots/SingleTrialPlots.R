# Code used to generate single trial plots

# Note: Since no seed was used, the plots will not match the plots used in the 
#       presentation exactly


source("../trial_design/base_trial_sims/3by3_function.R")
source("../trial_design/base_trial_sims/mtpi_function.R")
source("../trial_design/base_trial_sims/mtpi2_function.R")

library(ggplot2)




##################Code to Plot trial results using an mTPI design###############

#Function Parameters can be changed to explore design behavior under different toxicity assumptions
mtpi <- mTPI_sim_func(0.3, 0.05, c(10, 20, 40, 80, 160, 320), c(0.01, 0.05, 0.10, 0.25, 0.45, 0.55), 30, 0.40)

mtpi_plot <- ggplot(data = mtpi$patient_data, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(DLT), shape = decision), size = 3)+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  scale_shape_manual(
    name = "Decision",
    values = c("escalate" = 17, "stay" = 16, "de-escalate" = 25, "stay (at max)" = 21),
    labels = c("escalate" = "escalate", "stay" = "stay", "de-escalate" = "de-escalate"))+
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))+
  theme_classic()+
  labs(title = "Simulated mTPI Trial", x = "Patient", y = "Dose (mg)")

################################################################################



#################Code to Plot trial results using an mTPI2 design###############

#Function Parameters can be changed to explore design behavior under different toxicity assumptions
mtpi2 <- mTPI2_sim_func(0.3, 0.05, 0.05, c(10, 20, 40, 80, 160, 320), c(0.01, 0.05, 0.10, 0.25, 0.45, 0.55), 30, 0.40)

mtpi2_plot <- ggplot(data = mtpi2$patient_data, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(DLT), shape = decision), size = 3)+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40))+
  geom_line(color = "black", linetype = "dashed", lwd = 0.5)+
  scale_color_manual(name = "Toxicity Outcome",
                     values = c("0" = "forestgreen", "1" = "darkred"),
                     labels = c("0" = "No DLT", "1" = "DLT"))+
  scale_shape_manual(
    name = "Decision",
    values = c("escalate" = 17, "stay" = 16, "de-escalate" = 25, "stay (at max)" = 21),
    labels = c("escalate" = "escalate", "stay" = "stay", "de-escalate" = "de-escalate"))+
  theme_classic()+
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))+
  labs(title = "Simulated mTPI-2 Trial", x = "Patient", y = "Dose (mg)")

################################################################################



#################Code to Plot trial results using a 3+3  design#################

#Function Parameters can be changed to explore design behavior under different toxicity assumptions
threethree <- ThreeThree_sim_func(c(10, 20, 40, 80, 160, 320), c(0.01, 0.05, 0.10, 0.25, 0.45, 0.55))

threethree_plot <- ggplot(data = threethree$ThreeThree_data, aes(x = patient, y = dose))+
  geom_point(aes(color = factor(dlt)), size = 3)+
  scale_x_continuous(breaks = seq(3, nrow(threethree$ThreeThree_data), by = 3))+
  scale_y_continuous(breaks = c(10, 20, 40, 80, 160, 320))+
  geom_line(color = "black", linetype = "dashed")+
  scale_color_manual(
    values = c("1" = "darkred", "0" = "forestgreen"),
    labels = c("1" = "DLT", "0" = "No DLT")
  )+
  theme_classic()+
  labs(
    title = "3+3 Design Simulation",
    x = "Patient",
    y = "Dose (mg)"
  )

################################################################################

mtpi_plot
mtpi2_plot
threethree_plot
