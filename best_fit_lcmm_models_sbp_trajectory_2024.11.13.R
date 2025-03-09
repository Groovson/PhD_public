library(lcmm)
library(dplyr)
library(ggplot2)
library(ggplot2)

# Ensure group3 is a factor with distinct levels for the classes
# Convert numeric values to factor with the desired levels and labels

### Model 2 ###################################################################
hpt_lcmm2$group2a <- factor(hpt_lcmm2$group2,
                            levels = c(1, 2),      # Original numeric values
                            labels = c("Class 1", "Class 2"))  # New labels


# Plot with a legend and distinct colors for each latent class
p2a <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group2), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group2a, colour = group2a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red"))  # Custom colors for the classes
p2a


ggsave("model2_traj_lcmm_2024.11.20.png", p2a, width=15, height=13, units="cm", dpi=200, bg='white')

### Model 3 ###################################################################
hpt_lcmm2$group3a <- factor(hpt_lcmm2$group3,
                           levels = c(1, 2, 3),      # Original numeric values
                           labels = c("Class 1", "Class 2", "Class 3"))  # New labels


# Plot with a legend and distinct colors for each latent class
p3a <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group3a, colour = group3a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red", 
                                 "Class 3" = "green"))  # Custom colors for the classes
p3a


ggsave("model3_traj_lcmm_2024.11.20.png", p3a, width=15, height=13, units="cm", dpi=200, bg='white')


### Model 4 ##################################################################
hpt_lcmm2$group4a <- factor(hpt_lcmm2$group4,
                            levels = c(1, 2, 3, 4),      # Original numeric values
                            labels = c("Class 1", "Class 2", "Class 3",  "Class 4"))  # New labels


# Plot with a legend and distinct colors for each latent class
p4 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group4), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group4a, colour = group4a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red",
                                 "Class 3" = "purple",
                                 "Class 4" = "green"))  # Custom colors for the classes
p4

ggsave("model4_traj_lcmm_linear_20250203.png", p4, width=15, height=13, units="cm", dpi=200, bg='white')


### Model 5 ##################################################################
hpt_lcmm2$group5a <- factor(hpt_lcmm2$group5,
                            levels = c(1, 2, 3, 4, 5),      # Original numeric values
                            labels = c("Class 1", "Class 2", "Class 3",  "Class 4", "Class 5"))  # New labels


# Plot with a legend and distinct colors for each latent class
p5a <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group5a, colour = group5a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red",
                                 "Class 3" = "green",
                                 "Class 4" = "purple",
                                 "Class 5" = "pink"))  # Custom colors for the classes
p5a

ggsave("model5_traj_lcmm_2024.11.20.png", p5a, width=15, height=13, units="cm", dpi=200, bg='white')


### Model 6 ##################################################################
hpt_lcmm2$group6a <- factor(hpt_lcmm2$group6,
                            levels = c(1, 2, 3, 4, 5, 6),      # Original numeric values
                            labels = c("Class 1", "Class 2", "Class 3",  "Class 4", "Class 5", "Class 6"))  # New labels


# Plot with a legend and distinct colors for each latent class
p6a <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group6), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group6a, colour = group6a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red",
                                 "Class 3" = "green",
                                 "Class 4" = "pink",
                                 "Class 5" = "purple",
                                 "Class 6" = "yellow"))  # Custom colors for the classes
p6a

ggsave("model6_traj_lcmm_2024.11.20.png", p6a, width=15, height=13, units="cm", dpi=200, bg='white')


### Model 7 ##################################################################
hpt_lcmm2$group7a <- factor(hpt_lcmm2$group6,
                            levels = c(1, 2, 3, 4, 5, 6, 7),      # Original numeric values
                            labels = c("Class 1", "Class 2", "Class 3",  "Class 4", "Class 5", "Class 6", "Class 7"))  # New labels


# Plot with a legend and distinct colors for each latent class
p7a <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group7), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group7a, colour = group7a), method = "gam", 
              formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme_minimal() +
  theme(legend.position = "right") + 
  scale_colour_manual(values = c("Class 1" = "blue", 
                                 "Class 2" = "red",
                                 "Class 3" = "green",
                                 "Class 4" = "pink",
                                 "Class 5" = "purple",
                                 "Class 6" = "yellow"))  # Custom colors for the classes
p7a

ggsave("model7_traj_lcmm_2024.11.20.png", p7a, width=15, height=13, units="cm", dpi=200, bg='white')




p3 <- ggplot(hpt_lcmm2, aes(year2, systolic_bp, group = n_eid_14631)) +
  geom_smooth(aes(group = n_eid_14631, colour = group3), linewidth = 0.5, se = FALSE) + 
  geom_smooth(aes(group = group3), method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 2, se = TRUE) +
  scale_y_continuous(limits = c(80, 240)) + 
  labs(x = "Years", y = "Systolic BP", colour = "Latent Class", title = "Smoothed") +
  theme(legend.position = "none")
p3