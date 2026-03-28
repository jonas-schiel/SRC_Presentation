
source(here("../trial_design/LargeTrialSim.R"))

test <- get_data(doses = c(10, 20, 40, 80, 160, 320), dlt_probs = c(0.01, 0.05, 0.10, 0.25, 0.45, 0.55), 0.3, 0.05, 0.05, 0.45, 30, 1000)
library(ggplot2)

p1 <- ggplot(test$threethree, aes(x=MTD))+
  geom_histogram(fill = "yellow2", col = "black")+
  labs(
    title = "mTPI MTD Selection",
    x = "Selected MTD",
    y = "Count"
  )+
  theme_classic()

p2 <- ggplot(test$mtpi, aes(x=MTD))+
  geom_histogram(fill = "red2", col = "black")+
  labs(
    title = "mTPI-2 MTD Selection",
    x = "Selected MTD",
    y = "Count"
  )+
  theme_classic()

p3 <- ggplot(test$mtpi2, aes(x=MTD))+
  geom_histogram(fill = "blue2", col = "black")+
  labs(
    title = "3+3 MTD Selection",
    x = "Selected MTD",
    y = "Count"
  )+
  theme_classic()


combined <- rbind(
  transform(test$threethree, method = "3+3"),
  transform(test$mtpi, method = "mTPI"),
  transform(test$mtpi2, method = "mTPI2")
)

hist_plot <- ggplot(combined, aes(x = factor(MTD), fill = method)) +
  geom_bar(position = "dodge", alpha = 0.7, col = "black") +
  scale_fill_manual(values = c("3+3" = "seagreen3", "mTPI" = "wheat", "mTPI2" = "steelblue2")) +
  labs(title = "Selected Doses Across Models", subtitle = "10,000 Simulations using same dosing properties", x = "MTD", y = "Count", fill = "Method") +
  theme_classic()

p1
p2
p3
hist_plot