

# Required libraries
library(ggplot2)

library(dplyr);library(cowplot); library(tidyr)


## ------------------- Plot section for 10 stage high censoring ---------------##


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
  Value = c(pt300_10stage_hi_observed, pt300_10stage_hi_IHsurvrf),
  n = rep("n = 300", 2*length(pt300_10stage_hi_observed)),
  design = rep("10 stages: High Censoring Rate", 2*length(pt300_10stage_hi_observed)),
  setting = rep("2 Strata", 2*length(pt300_10stage_hi_observed))
)

## Plot using ggplot2 with ticks on Y-axis and labels for mean points
#highcens_300pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 300, 200 simulations, 10 stages, high censoring, 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(520, 620, by = 20)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt300_10stage_hi_observed), mean(pt300_10stage_hi_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6) +   # Label mean points
#  theme(legend.position = "none")


highcens_300pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(520, 590, by = 10), limits = c(520, 590)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt300_10stage_hi_observed), mean(pt300_10stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) + 
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt300_10stage_hi_observed), mean(pt300_10stage_hi_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 




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
  Value = c(pt500_10stage_hi_observed, pt500_10stage_hi_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_10stage_hi_observed)),
  design = rep("10 stages: High Censoring Rate", 2*length(pt500_10stage_hi_observed)),
  setting = rep("2 Strata", 2*length(pt500_10stage_hi_observed))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#highcens_500pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 500, 200 simulations, 10 stages, high censoring, 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(520, 620, by = 20)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt500_10stage_hi_observed), mean(pt500_10stage_hi_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6) +   # Label mean points
#  theme(legend.position = "none")

highcens_500pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  scale_y_continuous(breaks = seq(520, 590, by = 10), limits = c(520, 590))  +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_10stage_hi_observed), mean(pt500_10stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  ## suppress printing of y-axis labels
  theme(axis.title.y = element_blank(),  # Removes the y-axis title
        axis.text.y = element_blank(),    # Removes the y-axis labels
        axis.ticks.y = element_blank()) +  # Keeps the y-axis ticks 
  
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_10stage_hi_observed), mean(pt500_10stage_hi_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 

##################################################################
## DONE: 1 strata: plots for 300 pts, 10 stages, high censoring ##
## observational,                               ##
##################################################################
a <- read.csv("~/survrf/Outputs/1STRATA_29Aug_10stage__300pt_highcens_V1")
b <- read.csv("~/survrf/Outputs/1STRATA_04sep_10stage_300pt_highcens_V2")
c <- read.csv("~/survrf/Outputs/1STRATA_10stage_highcens_V3")
d <- read.csv("~/survrf/Outputs/1STRATA_10stage_300pt_highcens_V4")

  


## only 52 sim done: we want to make sure there are 200 data points in total, so if there aren't enough, we just run extra sims
## only 168 sims done in V2
pt300_10stage_hi_1_strat_observed <- c(a$observed, b$observed, c$observed, d$observed)[-1]

pt300_10stage_hi_1strat_IHsurvrf <- c(a$IHsurvrf, b$IHsurvrf, c$observed, d$observed)[-1]


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt300_10stage_hi_1_strat_observed)), rep("IHsurvrf", length(pt300_10stage_hi_1strat_IHsurvrf))),
  Value = c(pt300_10stage_hi_1_strat_observed, pt300_10stage_hi_1strat_IHsurvrf),
  n = rep("n = 300", 2*length(pt300_10stage_hi_1_strat_observed)),
  design = rep("10 stages: High Censoring Rate", 2*length(pt300_10stage_hi_1_strat_observed)),
  setting = rep("1 Strata", 2*length(pt300_10stage_hi_1_strat_observed))
)

highcens_300pt_10stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  scale_y_continuous(breaks = seq(520, 590, by = 10), limits = c(520, 590))  +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt300_10stage_hi_1_strat_observed), mean(pt300_10stage_hi_1strat_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt300_10stage_hi_1_strat_observed), mean(pt300_10stage_hi_1strat_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 

##################################################################
## DONE: 1 strata: plots for 500 pts, 10 stages, high censoring ##
## observational,                               ##
##################################################################
a <- read.csv("~/survrf/Outputs/1STRATA_01Sep_10stage_500pt_hicens_all")
b <- read.delim2("~/survrf/Outputs/done_running/1STRATA_500pt_hicens_obs_V2.txt", header=FALSE, comment.char="#")


# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V2 <- which(apply(b, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V2 <- data.frame()

# Loop through each matching index

for (index in matching_indices_V2) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V2 <- rbind(selected_rows_V2, b[index:min(index+1, nrow(b)), ])
}

## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
split_values_V2 <- strsplit(as.character(selected_rows_V2$X.1.204.544.9784.564.1007.555.5596.577.8582......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V2 <- lapply(split_values_V2, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[3])

# Convert to numeric if the values are numeric
observed2 <- as.numeric(third_values_V2)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V2 <- sapply(cleaned_split_values_V2, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf2 <- as.numeric(fourth_values_V2)


## only 52 sim done: we want to make sure there are 200 data points in total, so if there aren't enough, we just run extra sims
## only 168 sims done in V2
pt500_10stage_hi_1_strat_observed <- c(a$observed, observed2)[-1]

pt500_10stage_hi_1strat_IHsurvrf <- c(a$IHsurvrf, IHsurvrf2)[-1]


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_10stage_hi_1_strat_observed)), rep("IHsurvrf", length(pt500_10stage_hi_1strat_IHsurvrf))),
  Value = c(pt500_10stage_hi_1_strat_observed, pt500_10stage_hi_1strat_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_10stage_hi_1_strat_observed)),
  design = rep("10 stages: High Censoring Rate", 2*length(pt500_10stage_hi_1_strat_observed)),
  setting = rep("1 Strata", 2*length(pt500_10stage_hi_1_strat_observed))
)

highcens_500pt_10stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  scale_y_continuous(breaks = seq(520, 590, by = 10), limits = c(520, 590))  +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_10stage_hi_1_strat_observed), mean(pt500_10stage_hi_1strat_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  # Adding the black dot for the mean
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_10stage_hi_1_strat_observed), mean(pt500_10stage_hi_1strat_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 


####### grid plot

# Combine plots with a common label
plot_grid(highcens_300pt_10stage_2strata + guides(fill = FALSE), 
                           highcens_500pt_10stage_2strata + guides(fill = FALSE), 
          highcens_300pt_10stage_1strata,
          highcens_500pt_10stage_1strata,
                           align = "h", nrow = 2, ncol = 2)




## ----------------------------------------------------------------------------##

## ------------------- Plot section for 10 stage low censoring ---------------##

##################################################################
## DONE: plots for 300 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/300pt_lowcens_obs_V1.txt", header=FALSE, comment.char="#")
b <- read.delim2("~/survrf/Outputs/300pt_lowcens_obs_V2.txt", header=FALSE, comment.char="#")
c <- read.csv("~/survrf/Outputs/2STRATA_29Aug_10stage_300pt_lowcens_V3")
  
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

## we had 201 sims so we do -1
pt300_10stage_low_observed <- c(observed1, observed2, c$observed)[-1]

pt300_10stage_low_IHsurvrf <- c(IHsurvrf1, IHsurvrf2, c$IHsurvrf)[-1]



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt300_10stage_low_observed)), rep("IHsurvrf", length(pt300_10stage_low_IHsurvrf))),
  Value = c(pt300_10stage_low_observed, pt300_10stage_low_IHsurvrf),
  n = rep("n = 300", 2*length(pt300_10stage_low_IHsurvrf)),
  design = rep("10 stages: Low Censoring Rate", 2*length(pt300_10stage_low_observed)),
  setting = rep("2 Strata", 2*length(pt300_10stage_low_observed))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#lowcens_300pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 300, 185 simulations, 10 stages, moderate censoring (30%), 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(440, 540, by = 10)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt300_10stage_low_observed), mean(pt300_10stage_low_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points
#

lowcens_300pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(450, 520, by = 10), limits = c(450, 520)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt300_10stage_low_observed), mean(pt300_10stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt300_10stage_low_observed), mean(pt300_10stage_low_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 




##################################################################
## DONE-- plots for 500 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/500pt_lowcens_run1.txt", header=FALSE, comment.char="#")
b <- read.delim2("~/survrf/Outputs/500pt_lowcens_run2.txt", header=FALSE, comment.char="#")
c <- read.csv("~/survrf/Outputs/2STRATA_31Aug_10stage_500pt_lowcens_run3")


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

pt500_10stage_low_observed <- c(observed1, observed2, c$observed)

pt500_10stage_low_IHsurvrf <- c(IHsurvrf1, IHsurvrf2, c$IHsurvrf)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_10stage_low_observed)), rep("IHsurvrf", length(pt500_10stage_low_IHsurvrf))),
  Value = c(pt500_10stage_low_observed, pt500_10stage_low_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_10stage_low_observed)),
  design = rep("10 stages: Low Censoring Rate", 2*length(pt500_10stage_low_observed)),
  setting = rep("2 Strata", 2*length(pt500_10stage_low_observed))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#lowcens_500pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 500, 153 simulations, 10 stages, moderate censoring (26%), 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(440, 540, by = 10)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt500_10stage_low_observed), mean(pt500_10stage_low_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points
#


lowcens_500pt_10stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(450, 520, by = 10), limits = c(450, 520)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_10stage_low_observed), mean(pt500_10stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_10stage_low_observed), mean(pt500_10stage_low_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



##################################################################
## DONE-- 1 strata: plots for 300 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

one_strat_10stage_lowcensa <- read.csv("~/survrf/Outputs/1STRATA_03sep_10stage_300pt_lowcens_V1")
one_strat_10stage_lowcensb <- read.csv("~/survrf/Outputs/1STRATA_01Sep_10stage_300pt_lowcens_V2")

one_strat_10stage_lowcens_OBS <- c(one_strat_10stage_lowcensa$observed, one_strat_10stage_lowcensb$observed)[-1]

one_strat_10stage_lowcens_IHsurvrf <- c(one_strat_10stage_lowcensa$IHsurvrf, one_strat_10stage_lowcensb$IHsurvrf)[-1]

# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(one_strat_10stage_lowcens_OBS)), rep("IHsurvrf", length(one_strat_10stage_lowcens_IHsurvrf))),
  Value = c(one_strat_10stage_lowcens_OBS, one_strat_10stage_lowcens_IHsurvrf),
  n = rep("n = 300", 2*length(one_strat_10stage_lowcens_IHsurvrf)),
  design = rep("10 stages: Low Censoring Rate", 2*length(one_strat_10stage_lowcens_IHsurvrf)),
  setting = rep("1 Strata", 2*length(one_strat_10stage_lowcens_IHsurvrf))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 300, ONE STRATA, 100 simulations, 10 stages, low censoring, 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(530, 580, by = 10)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(one_strat_10stage_lowcens_OBS), mean(one_strat_10stage_lowcens_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points
#

lowcens_300pt_10stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(450, 520, by = 10), limits = c(450, 520)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(one_strat_10stage_lowcens_OBS), mean(one_strat_10stage_lowcens_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(one_strat_10stage_lowcens_OBS), mean(one_strat_10stage_lowcens_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



##################################################################
## 1 strata: plots for 500 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

one_strat_10stage_lowcens_500a <- read.csv("~/survrf/Outputs/1STRATA_29Aug_10stage_500pt_lowcens_V1")
one_strat_10stage_lowcens_500b <- read.csv("~/survrf/Outputs/1STRATA_500pt_locens_201-220")


######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/done_running/1STRATA_500pt_locens_obs_V2.txt", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}


## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1.109.480.5316.505.3315.491.6391.503.0702......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)

one_strat_10stage_lowcens_500_OBS <- c(one_strat_10stage_lowcens_500a$observed, observed1, one_strat_10stage_lowcens_500b$observed)[-1]

one_strat_10stage_lowcens_500_IHsurvrf <- c(one_strat_10stage_lowcens_500a$IHsurvrf, IHsurvrf1, one_strat_10stage_lowcens_500b$IHsurvrf)[-1]

# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(one_strat_10stage_lowcens_500_OBS)), rep("IHsurvrf", length(one_strat_10stage_lowcens_500_IHsurvrf))),
  Value = c(one_strat_10stage_lowcens_500_OBS, one_strat_10stage_lowcens_500_IHsurvrf),
  n = rep("n = 500", 2*length(one_strat_10stage_lowcens_500_IHsurvrf)),
  design = rep("10 stages: Low Censoring Rate", 2*length(one_strat_10stage_lowcens_500_IHsurvrf)),
  setting = rep("1 Strata", 2*length(one_strat_10stage_lowcens_500_IHsurvrf))
)


lowcens_500pt_10stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(450, 520, by = 10), limits = c(450, 520)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(one_strat_10stage_lowcens_500_OBS), mean(one_strat_10stage_lowcens_500_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(one_strat_10stage_lowcens_500_OBS), mean(one_strat_10stage_lowcens_500_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



####### grid plot

# Combine plots with a common label
plot_grid(lowcens_300pt_10stage_2strata + guides(fill = FALSE), 
          lowcens_500pt_10stage_2strata + guides(fill = FALSE), 
          lowcens_300pt_10stage_1strata,
          lowcens_500pt_10stage_1strata,
          align = "h", nrow = 2, ncol = 2)

## ----------------------------------------------------------------------------##

## ------------------- Plot section for 25 stage  ---------------##

##################################################################
## plots for 500 pts, 25 stages, low censoring ##
## observational,                               ##
##################################################################

## V1:only 7 sim done
## V2: only 107 sim done
## V3: only 12 sim done
## V4: only 22 sim done
## V5: only 30 sim done
## V6: only 51 sim done

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/done_running/500pt_25stage_locens_obs_V7.txt", header=FALSE, comment.char="#")
b <- read.delim2("~/survrf/Outputs/done_running/500pt_25stage_lowcens_376:531.txt", header=FALSE, comment.char="#")


# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))
matching_indices_V2 <- which(apply(b, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}

selected_rows_V2 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V2) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V2 <- rbind(selected_rows_V2, b[index:min(index+1, nrow(b)), ])
}


## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1.205..1137.04.1315.774.1315.179.1303.026......NA......................0.), " ")
split_values_V2 <- strsplit(as.character(selected_rows_V2$X.1.376.1163.111.1196.877.1290.115.1316.58......NA......................0.), " ")

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

pt500_25stage_low_observed <- c(1179.303, 1174.092, 1186.392, 1168.456, 1164.47,
                                1173.317, 1159.486, 1135.751, 1172.255, 1170.847,
                                1165.927, 1180.468, 1161.419, 1167.43, 1188.548, 1163.283,
                                1166.918, 1167.967, 1165.285, 1153.544, 1158.596,
                                1159.414, 1172.082, 1150.766, 1162.323, observed1, observed2)

pt500_25stage_low_IHsurvrf <- c(1228.147, 1157.587, 1228.418, 1293.653, 1167.892,
                                1160.929, 1197.294, 1240.757, 1254.862, 1217.989,
                                1296.912, 1297.104, 1301.163, 1304.961, 1252.483, 1329.84,
                                1241.718, 1253.786, 1306.99, 1320.823, 1289.146,
                                1198.282, 1238.886, 1270.393, 1260.407, IHsurvrf1, IHsurvrf2)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25stage_low_observed)), rep("IHsurvrf", length(pt500_25stage_low_IHsurvrf))),
  Value = c(pt500_25stage_low_observed, pt500_25stage_low_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_25stage_low_IHsurvrf)),
  design = rep("25 stages: Low Censoring Rate", 2*length(pt500_25stage_low_IHsurvrf)),
  setting = rep("2 Strata", 2*length(pt500_25stage_low_IHsurvrf))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#lowcens_500pt_25stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 500, 15 simulations, 25 stages, moderate censoring (37%), 
#       Observational setting, 10,000 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(1100, 1300, by = 20)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt500_25stage_low_observed), mean(pt500_25stage_low_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6)   # Label mean points

lowcens_500pt_25stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(1120, 1350, by = 30), limits = c(1120, 1350)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25stage_low_observed), mean(pt500_25stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_25stage_low_observed), mean(pt500_25stage_low_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



##################################################################
## plots for 500 pts, 25 stages, high censoring ##
## observational,                               ##
##################################################################

two_strat_25stage_hicens_1 <- read.csv("~/survrf/Outputs/2STRATA_01Sep_25stage_500pt_hicens_PAR1")

######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/done_running/500pt_25stage_hicens_201:372.txt", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}


## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1.205.1718.032.1831.851.1886.584.1857.042......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)

## using a threshold of 0.25
## 5 sims done
## 105 sims done
pt500_25stage_hi_observed <- c(1707.919, 1691.335, 1675.723, 1703.018,
                               1695.49, 1691.993, 1692.666, 1713.72,
                               1711.027, 1712.989, two_strat_25stage_hicens_1$observed, observed1)
pt500_25stage_hi_IHsurvrf <- c(1808.825, 1817.265, 1747.049, 1723.135,
                               1790.413, 1816.506, 1844.138, 1757.585,
                               1860.377, 1806.534, two_strat_25stage_hicens_1$IHsurvrf, IHsurvrf1)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25stage_hi_observed)), rep("IHsurvrf", length(pt500_25stage_hi_IHsurvrf))),
  Value = c(pt500_25stage_hi_observed, pt500_25stage_hi_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_25stage_hi_IHsurvrf)),
  design = rep("25 stages: High Censoring Rate", 2*length(pt500_25stage_hi_IHsurvrf)),
  setting = rep("2 Strata", 2*length(pt500_25stage_hi_IHsurvrf))
)

# Plot using ggplot2 with ticks on Y-axis and labels for mean points
#highcens_500pt_25stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
#  geom_boxplot(alpha = 0, outlier.shape = NA) +  # Make quartiles more visible
#  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.5) +  # Plot points with same color as boxplot and make them more transparent
#  geom_point(stat = "summary", fun = "mean", shape = 19, size = 3, color = "black") +  # Plot mean points
#  labs(title = "Sample size of 500, 10 simulations, 25 stages, higher censoring (50%), 
#       Observational setting, 100 sample for value evaluation",
#       x = NULL, y = "Values") +
#  theme_minimal() +
#  scale_color_manual(values = c("observed" = "blue", "IHsurvrf" = "red")) +  # Match point colors to boxplot colors
#  scale_y_continuous(breaks = seq(1650, 1850, by = 20)) +  # Set ticks on Y-axis by 10
#  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
#                              Value = c(mean(pt500_25stage_hi_observed), mean(pt500_25stage_hi_IHsurvrf))),
#            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = -0.5,color = "black", size = 6) +  # Label mean points
#  
#  theme(legend.position = "none")

highcens_500pt_25stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(1650, 1900, by = 30), limits = c(1650, 1900)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25stage_hi_observed), mean(pt500_25stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_25stage_hi_observed), mean(pt500_25stage_hi_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 


##################################################################
## 1 strata: plots for 500 pts, 25 stages, low censoring ##
## observational,                               ##
##################################################################


######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/done_runnning/1STRATA_25stage_500pt_locens_V1.txt", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}


## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1..3.1186.392.1210.587.1311.914.1287.359......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)

pt500_25_lo1_observed <- observed1
pt500_25_lo1_IHsurvrf <- IHsurvrf1


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25_lo1_observed)), rep("IHsurvrf", length(pt500_25_lo1_IHsurvrf))),
  Value = c(pt500_25_lo1_observed, pt500_25_lo1_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_25_lo1_observed)),
  design = rep("25 stages: Low Censoring Rate", 2*length(pt500_25_lo1_observed)),
  setting = rep("1 Strata", 2*length(pt500_25_lo1_observed))
)


locens_500pt_15stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(1150, 1300, by = 30), limits = c(1150, 1300)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25_lo1_observed), mean(pt500_25_lo1_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_25_lo1_observed), mean(pt500_25_lo1_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 
  



##################################################################
## 1 strata: plots for 500 pts, 25 stages, high censoring ##
## observational,                               ##
##################################################################


######### data preprocessing for the .txt files that didn't finish running
a <- read.delim2("~/survrf/Outputs/1STRATA_25stage_500pt_hicens_V1", header=FALSE, comment.char="#")

# Find row indices where the pattern "Estimation-True Optimal" appears in any column
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

# Initialize an empty vector to store the selected rows
selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}


## splitting the values of the second column which contain the numeric results of the sims
#########
######### NOTE: this subsetting the column name part most likely to make errors
#########
split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1..4.1691.993.1802.992.1861.475.1853.67......NA......................0.), " ")


# Remove any blank elements ("") from each split string
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])

## now, we subset the observed and the IHsurvrf

##### observed #######
## observed is column 2's 3rd number
# Extract the third value from each split string
third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])

# Convert to numeric if the values are numeric
observed1 <- as.numeric(third_values_V1)

##### IHsurvrf #######
## IHsurvrf is column 2's 4th number
fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])

# Convert to numeric if the values are numeric
IHsurvrf1 <- as.numeric(fourth_values_V1)

pt500_25_hi1_observed <- observed1
pt500_25_hi1_IHsurvrf <- IHsurvrf1


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_25_hi1_observed)), rep("IHsurvrf", length(pt500_25_hi1_IHsurvrf))),
  Value = c(pt500_25_hi1_observed, pt500_25_hi1_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_25_hi1_observed)),
  design = rep("25 stages: High Censoring Rate", 2*length(pt500_25_hi1_observed)),
  setting = rep("1 Strata", 2*length(pt500_25_hi1_observed))
)


hicens_500pt_15stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(1650, 1890, by = 30), limits = c(1650, 1890)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_25_hi1_observed), mean(pt500_25_hi1_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_25_hi1_observed), mean(pt500_25_hi1_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



####### grid plot

# Combine plots with a common label
plot_grid(lowcens_500pt_25stage_2strata + guides(fill = FALSE), 
          highcens_500pt_25stage_2strata + guides(fill = FALSE),
          hicens_500pt_15stage_1strata, 
          locens_500pt_15stage_1strata,
          align = "h", nrow = 2, ncol = 2)

## ---------------------------------------------------------------##










