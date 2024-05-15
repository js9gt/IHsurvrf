

# Required libraries
library(ggplot2)



##################################################################
## plots for 75 patients, 20 sims, 30/40 stages, ~10% censoring ##
## observational, 100 samples for value estimation.             ##
##################################################################

pt75observed <- c(302.2944, 360.6497, 362.3956, 293.4051, 386.3928, 372.6228, 311.3882, 344.9238, 352.8118, 346.8968,
                  330.3893, 357.0209, 362.0792, 337.8015, 347.5771, 324.2894, 347.6283, 346.7461, 334.6861, 348.1616)

pt75IHsurvrf <- c(283.7597, 354.7568, 348.1319, 307.8553, 348.8683, 306.6561, 304.4959, 364.1864, 286.1369, 377.5856,
                  319.161, 354.7152, 373.7673, 333.3085, 264.6007, 358.9236, 292.7851, 372.0944, 271.1077, 263.4456)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt75observed)), rep("IHsurvrf", length(pt75IHsurvrf))),
  Value = c(pt75observed, pt75IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 75, 20 simulations, 30/40 stages, ~10% censoring, 
       Observational setting, 100 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(250, 400, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt75observed), mean(pt75IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points


###################################################################
## plots for 200 patients, 20 sims, 30/40 stages, ~10% censoring ##
## observational, 100 samples for value estimation.              ##
###################################################################


pt200observed <-c(302.2944, 360.6497, 362.3956, 293.4051, 386.3928, 372.6228, 311.3882, 344.9238, 352.8118, 346.8968,
                  330.3893, 357.0209, 362.0792, 337.8015, 347.5771, 324.2894, 347.6283, 346.7461, 334.6861, 348.1616)
pt200IHsurvrf <- c(326.4924, 312.2935, 386.7455, 316.7066, 347.3867, 304.5746, 267.7904, 386.0886, 348.1822, 356.3409,
                   280.0495, 376.1314, 263.2632, 270.3984, 317.2554, 334.897, 315.9025, 334.5974, 349.057, 343.8445)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt200observed)), rep("IHsurvrf", length(pt200IHsurvrf))),
  Value = c(pt200observed, pt200IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 200, 20 simulations, 30/40 stages, ~10% censoring, 
       Observational setting, 100 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(250, 400, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt200observed), mean(pt200IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points



