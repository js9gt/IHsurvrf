

# Required libraries
library(ggplot2)



##################################################################
## plots for 500 pts, 25 stages, high censoring ##
## observational,                               ##
##################################################################

## using a threshold of 0.3
## 2 sims done
## 102 sims done
pt500_25stage_hi_observed <- c(1707.919, 1691.335, 1675.723, 1703.018)
pt500_25stage_hi_IHsurvrf <- c(1808.825, 1817.265, 1747.049, 1723.135)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25stage_hi_observed)), rep("IHsurvrf", length(pt500_25stage_hi_IHsurvrf))),
  Value = c(pt500_25stage_hi_observed, pt500_25stage_hi_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 500, 4 simulations, 25 stages, higher censoring (50%), 
       Observational setting, 100 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(1650, 1850, by = 20)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25stage_hi_observed), mean(pt500_25stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points
  


##################################################################
## plots for 500 pts, 25 stages, low censoring ##
## observational,                               ##
##################################################################

## V1:only 7 sim done
## V2: only 107 sim done
## V3: only 1 sim done
## V4
pt500_25stage_low_observed <- c(1179.303, 1174.092, 1186.392, 1168.456, 1164.47,
                                1173.317, 1159.486, 1135.751, 1172.255, 1170.847,
                                1165.927, 1180.468, 1161.419, 1167.43, 1188.548)

pt500_25stage_low_IHsurvrf <- c(1228.147, 1157.587, 1228.418, 1293.653, 1167.892,
                                1160.929, 1197.294, 1240.757, 1254.862, 1217.989,
                                1296.912, 1297.104, 1301.163, 1304.961, 1252.483)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25stage_low_observed)), rep("IHsurvrf", length(pt500_25stage_low_IHsurvrf))),
  Value = c(pt500_25stage_low_observed, pt500_25stage_low_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 500, 15 simulations, 25 stages, moderate censoring (37%), 
       Observational setting, 10,000 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(1100, 1300, by = 20)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25stage_low_observed), mean(pt500_25stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points


##################################################################
## plots for 300 pts, 10 stages, high censoring ##
## observational--- DONE---                     ##
##################################################################

a <- read.csv("~/survrf/Outputs/08Aug_300pt_hicens_obs_RESUTS_V1")
b <- read.csv("~/survrf/Outputs/08Aug_300pt_hicens_obs_RESUTS_V2")
comb <- rbind(a[85:185, ], b[186:286, ])

## only 58 sim done: we want to make sure there are 200 data points in total, so if there aren't enough, we just run extra sims
## run1: 86 sim done
## run2: 183 sims done
pt300_10stage_hi_observed <- comb$observed
pt300_10stage_hi_IHsurvrf <- comb$IHsurvrf


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt300_10stage_hi_observed)), rep("IHsurvrf", length(pt300_10stage_hi_IHsurvrf))),
  Value = c(pt300_10stage_hi_observed, pt300_10stage_hi_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 300, 200 simulations, 10 stages, high censoring, 
       Observational setting, 10,000 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(520, 600, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt300_10stage_hi_observed), mean(pt300_10stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points


##################################################################
## plots for 300 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/300pt_lowcens_obs_V1.txt", header=FALSE, comment.char="#")
b <- read.delim2("~/survrf/Outputs/300pt_lowcens_obs_V2.txt", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))
matching_indices_V2 <- which(apply(b, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()
selected_rows_V2 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}

for (index in matching_indices_V2) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V2 <- rbind(selected_rows_V2, b[index:min(index+1, nrow(b)), ])
}

## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1..1.475.5547..483.679.487.9192.511.5064......NA......................0.), " ")
split_values_V2 <- strsplit(as.character(selected_rows_V2$X.101.NA.473.1292.491.7479.487.7084.504.2211......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])
cleaned_split_values_V2 <- lapply(split_values_V2, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])
third_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)
observed2 <- as.numeric(third_values_V2)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])
fourth_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)
IHsurvrf2 <- as.numeric(fourth_values_V2)

## only 79 sim done: we want to make sure there are 200 data points in total, so if there aren't enough, we just run extra sims
pt300_10stage_low_observed <- c(observed1, observed2)

pt300_10stage_low_IHsurvrf <- c(IHsurvrf1, IHsurvrf2)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt300_10stage_low_observed)), rep("IHsurvrf", length(pt300_10stage_low_IHsurvrf))),
  Value = c(pt300_10stage_low_observed, pt300_10stage_low_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 300, 89 simulations, 10 stages, moderate censoring (30%), 
       Observational setting, 10,000 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(440, 540, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt300_10stage_low_observed), mean(pt300_10stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points


##################################################################
## plots for 500 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/500pt_lowcens_run1.txt", header=FALSE, comment.char="#")
b <- read.delim2("~/survrf/Outputs/500pt_lowcens_run2.txt", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))
matching_indices_V2 <- which(apply(b, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()
selected_rows_V2 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}

for (index in matching_indices_V2) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V2 <- rbind(selected_rows_V2, b[index:min(index+1, nrow(b)), ])
}

## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1..1..474.618..469.668.489.2657.508.2545......NA......................0.), " ")
split_values_V2 <- strsplit(as.character(selected_rows_V2$X.101.NA.473.5441.494.5405.485.2962.509.0162......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])
cleaned_split_values_V2 <- lapply(split_values_V2, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])
third_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)
observed2 <- as.numeric(third_values_V2)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])
fourth_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)
IHsurvrf2 <- as.numeric(fourth_values_V2)

pt500_10stage_low_observed <- c(observed1, observed2)

pt500_10stage_low_IHsurvrf <- c(IHsurvrf1, IHsurvrf2)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_10stage_low_observed)), rep("IHsurvrf", length(pt500_10stage_low_IHsurvrf))),
  Value = c(pt500_10stage_low_observed, pt500_10stage_low_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 500, 67 simulations, 10 stages, moderate censoring (26%), 
       Observational setting, 10,000 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(440, 540, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_10stage_low_observed), mean(pt500_10stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points



##################################################################
## DONE: plots for 500 pts, 10 stages, high censoring ##
## observational,                               ##
##################################################################

a <- read.csv("~/survrf/Outputs/08Aug_500pt_hicens_obs_RESUTS_V2")
b <- read.csv("~/survrf/Outputs/08Aug_500pt_hicens_obs_RESUTS")
comb <- rbind(b[53:153, ], a[154:254, ])

## only 52 sim done: we want to make sure there are 200 data points in total, so if there aren't enough, we just run extra sims
## only 168 sims done in V2
pt500_10stage_hi_observed <- comb$observed

pt500_10stage_hi_IHsurvrf <- comb$IHsurvrf


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_10stage_hi_observed)), rep("IHsurvrf", length(pt500_10stage_hi_IHsurvrf))),
  Value = c(pt500_10stage_hi_observed, pt500_10stage_hi_IHsurvrf)
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
  labs(title = "Sample size of 500, 200 simulations, 10 stages, high censoring, 
       Observational setting, 10,000 sample for value evaluation",
       x = NULL, y = "Values") +
  theme_minimal() +
  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
  scale_y_continuous(breaks = seq(530, 580, by = 10)) +  # Set ticks on Y-axis by 10
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_10stage_hi_observed), mean(pt500_10stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points


