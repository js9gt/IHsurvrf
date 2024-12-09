

# Required libraries
library(ggplot2)

library(dplyr);library(cowplot); library(tidyr)


## ------------------- Plot section for 10 stage high censoring ---------------##


##################################################################
## plots for 300 pts, 10 stages, high censoring ##
## observational--- DONE---                     ##
##################################################################

two_strat_10stage_hicens1 <- read.csv("~/survrf/Outputs/10stage_300pt_hicens_1:50")
two_strat_10stage_hicens2 <- read.csv("~/survrf/Outputs/10stage_300pt_hicens_51:100")
two_strat_10stage_hicens3 <- read.csv("~/survrf/Outputs/10stage_300pt_hicens_101:100")
two_strat_10stage_hicens4 <- read.csv("~/survrf/Outputs/10stage_300pt_hicens_151:200")

pt300_10stage_hi_observed <- c(two_strat_10stage_hicens1$observed, two_strat_10stage_hicens2$observed,
                               two_strat_10stage_hicens3$observed, two_strat_10stage_hicens4$observed)
pt300_10stage_hi_IHsurvrf <- c(two_strat_10stage_hicens1$IHsurvrf, two_strat_10stage_hicens2$IHsurvrf,
                               two_strat_10stage_hicens3$IHsurvrf, two_strat_10stage_hicens4$IHsurvrf)

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


two_strat_10stage_500pt_hicens1 <- read.csv("~/survrf/Outputs/10stage_500pt_hicens_151:200")
two_strat_10stage_500pt_hicens2 <- read.csv("~/survrf/Outputs/10stage_500pt_hicens_51:100")
two_strat_10stage_500pt_hicens3 <- read.csv("~/survrf/Outputs/10stage_500pt_hicens_1:50")
two_strat_10stage_500pt_hicens4 <- read.csv("~/survrf/Outputs/10stage_500pt_hicens_101:150")


pt500_10stage_hi_observed <- c(two_strat_10stage_500pt_hicens1$observed, two_strat_10stage_500pt_hicens2$observed,
                               two_strat_10stage_500pt_hicens3$observed, two_strat_10stage_500pt_hicens4$observed)

pt500_10stage_hi_IHsurvrf <- c(two_strat_10stage_500pt_hicens1$IHsurvrf, two_strat_10stage_500pt_hicens2$IHsurvrf,
                               two_strat_10stage_500pt_hicens3$IHsurvrf, two_strat_10stage_500pt_hicens4$IHsurvrf)


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
onestrat_300pt_hicens1 <- read.csv("~/survrf/Outputs/1strata_10stage_300pt_hicens_101:150")
onestrat_300pt_hicens2 <- read.csv("~/survrf/Outputs/1strata_10stage_300pt_hicens_1:50")
onestrat_300pt_hicens3 <- read.csv("~/survrf/Outputs/1strata_10stage_300pt_hicens_51:100")
onestrat_300pt_hicens4 <- read.csv("~/survrf/Outputs/1strata_10stage_300pt_hicens_151:200")

  
pt300_10stage_hi_1_strat_observed <- c(onestrat_300pt_hicens1$observed, onestrat_300pt_hicens2$observed,
                                       onestrat_300pt_hicens3$observed, onestrat_300pt_hicens4$observed)

pt300_10stage_hi_1strat_IHsurvrf <- c(onestrat_300pt_hicens1$IHsurvrf, onestrat_300pt_hicens2$IHsurvrf,
                                      onestrat_300pt_hicens3$IHsurvrf, onestrat_300pt_hicens4$IHsurvrf)


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
onestrat_500pt_hicens1 <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_hicens_1:50")
onestrat_500pt_hicens2 <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_hicens_101:150")
onestrat_500pt_hicens3 <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_hicens_51:100")
onestrat_500pt_hicens4 <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_hicens_151:200")
onestrat_500pt_hicens5 <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_hicens_201:207")



pt500_10stage_hi_1_strat_observed <- c(onestrat_500pt_hicens1$observed, onestrat_500pt_hicens2$observed,
                                       onestrat_500pt_hicens3$observed, onestrat_500pt_hicens4$observed,
                                       onestrat_500pt_hicens5$observed)

pt500_10stage_hi_1strat_IHsurvrf <- c(onestrat_500pt_hicens1$IHsurvrf, onestrat_500pt_hicens2$IHsurvrf,
                                      onestrat_500pt_hicens3$IHsurvrf, onestrat_500pt_hicens4$IHsurvrf,
                                      onestrat_500pt_hicens5$IHsurvrf)


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

# Save the plot with a higher resolution
ggsave("10stage_hicens.png", plot = last_plot(), dpi = 300, width = 10, height = 8)



### now looking at the standard deviations
## 2 strata: 
# ------ 300 high censoring
sd(pt300_10stage_hi_IHsurvrf)

# ------ 500 high censoring (smaller sd than 300 pt)
sd(pt500_10stage_hi_IHsurvrf)

## 1 strata:
# ------ 300 high censoring (smaller sd than 2 strata & higher performance)
sd(pt300_10stage_hi_1strat_IHsurvrf)

# ------ 500 high censoring (smaller sd than 2 strata 500 pt, but larger sd than 300)
sd(pt500_10stage_hi_1strat_IHsurvrf)

## ----------------------------------------------------------------------------##

## ------------------- Plot section for 10 stage low censoring ---------------##

##################################################################
## DONE: plots for 300 pts, 10 stages, low censoring ##
## observational,                               ##
##################################################################

tenstage_300pt_locens1 <- read.csv("~/survrf/Outputs/10stage_300pt_locens_1:50")
tenstage_300pt_locens2 <- read.csv("~/survrf/Outputs/10stage_300pt_locens_51:100")
tenstage_300pt_locens3 <- read.csv("~/survrf/Outputs/10stage_300pt_locens_101:150")
tenstage_300pt_locens4 <- read.csv("~/survrf/Outputs/10stage_300pt_locens_151:200")



## we had 201 sims so we do -1
pt300_10stage_low_observed <- c(tenstage_300pt_locens1$observed, tenstage_300pt_locens2$observed, 
                                tenstage_300pt_locens3$observed, tenstage_300pt_locens4$observed)

pt300_10stage_low_IHsurvrf <- c(tenstage_300pt_locens1$IHsurvrf, tenstage_300pt_locens2$IHsurvrf,
                                tenstage_300pt_locens3$IHsurvrf, tenstage_300pt_locens4$IHsurvrf)



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
  scale_y_continuous(breaks = seq(450, 530, by = 10), limits = c(450, 530)) +
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


lowcens_10stage_500pt_1 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_1:50")
lowcens_10stage_500pt_2 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_151:200")
lowcens_10stage_500pt_3 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_51:100")
lowcens_10stage_500pt_4 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_101:150")
lowcens_10stage_500pt_5 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_201:233")
lowcens_10stage_500pt_6 <- read.csv("~/survrf/Outputs/10stage_500pt_locens_234:243")

pt500_10stage_low_observed <- c(lowcens_10stage_500pt_1$observed, lowcens_10stage_500pt_2$observed,
                                lowcens_10stage_500pt_3$observed, lowcens_10stage_500pt_4$observed,
                                lowcens_10stage_500pt_5$observed, lowcens_10stage_500pt_6$observed)
pt500_10stage_low_IHsurvrf <- c(lowcens_10stage_500pt_1$IHsurvrf, lowcens_10stage_500pt_2$IHsurvrf,
                                lowcens_10stage_500pt_3$IHsurvrf, lowcens_10stage_500pt_4$IHsurvrf,
                                lowcens_10stage_500pt_5$IHsurvrf, lowcens_10stage_500pt_6$IHsurvrf)


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
  scale_y_continuous(breaks = seq(450, 530, by = 10), limits = c(450, 530)) +
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

one_strat_10stage_lowcensa <- read.csv("~/survrf/Outputs/1STRATA_10stage_300pt_locens_101:150")
one_strat_10stage_lowcensb <- read.csv("~/survrf/Outputs/1STRATA_10stage_300pt_locens_151:200")
one_strat_10stage_lowcensc <- read.csv("~/survrf/Outputs/1STRATA_10stage_300pt_locens_1:50")
one_strat_10stage_lowcensd <- read.csv("~/survrf/Outputs/1STRATA_10stage_300pt_locens_51:100")

one_strat_10stage_lowcens_OBS <- c(one_strat_10stage_lowcensa$observed, one_strat_10stage_lowcensb$observed,
                                   one_strat_10stage_lowcensc$observed, one_strat_10stage_lowcensd$observed)

one_strat_10stage_lowcens_IHsurvrf <- c(one_strat_10stage_lowcensa$IHsurvrf, one_strat_10stage_lowcensb$IHsurvrf,
                                        one_strat_10stage_lowcensc$IHsurvrf, one_strat_10stage_lowcensd$IHsurvrf)

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
  scale_y_continuous(breaks = seq(450, 530, by = 10), limits = c(450, 530)) +
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

one_strat_10stage_lowcens_500a <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_151:200")
one_strat_10stage_lowcens_500b <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_101:150")
one_strat_10stage_lowcens_500c <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_1:50")
one_strat_10stage_lowcens_500d <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_51:100")
one_strat_10stage_lowcens_500e <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_201:244")
one_strat_10stage_lowcens_500f <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_245:249")
one_strat_10stage_lowcens_500g <- read.csv("~/survrf/Outputs/1STRATA_10stage_500pt_locens_250")




one_strat_10stage_lowcens_500_OBS <- c(one_strat_10stage_lowcens_500a$observed, one_strat_10stage_lowcens_500b$observed,
                                       one_strat_10stage_lowcens_500c$observed, one_strat_10stage_lowcens_500d$observed,
                                       one_strat_10stage_lowcens_500e$observed, one_strat_10stage_lowcens_500f$observed)
one_strat_10stage_lowcens_500_IHsurvrf <- c(one_strat_10stage_lowcens_500a$IHsurvrf, one_strat_10stage_lowcens_500b$IHsurvrf,
                                            one_strat_10stage_lowcens_500c$IHsurvrf, one_strat_10stage_lowcens_500d$IHsurvrf,
                                            one_strat_10stage_lowcens_500e$IHsurvrf, one_strat_10stage_lowcens_500f$IHsurvrf)

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
  scale_y_continuous(breaks = seq(450, 530, by = 10), limits = c(450, 530)) +
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

# Save the plot with a higher resolution
ggsave("10stage_locens.png", plot = last_plot(), dpi = 300, width = 10, height = 8)

### now looking at the standard deviations
## 2 strata: 
# ------ 300 low censoring
sd(pt300_10stage_low_IHsurvrf)

# ------ 500 low censoring (larger sd than 300 pt)
sd(pt500_10stage_low_IHsurvrf)

## 1 strata:
# ------ 300 low censoring (smaller sd than 2 strata & higher performance)
sd(one_strat_10stage_lowcens_IHsurvrf)

# ------ 500 low censoring (smaller sd than 2 strata 500 pt, and smaller sd than the 1 strata 300 pt)
sd(one_strat_10stage_lowcens_500_IHsurvrf)

## ----------------------------------------------------------------------------##

## ------------------- Plot section for 15 stage  ---------------##

##################################################################
## plots for 500 pts, 15 stages, low censoring ##
## observational,                               ##
##################################################################


two_strat_15stage_locens_1 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_26:50")
two_strat_15stage_locens_2 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_76:100")
two_strat_15stage_locens_3 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_101:125")
two_strat_15stage_locens_4 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_1:25")
two_strat_15stage_locens_5 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_151:175")
two_strat_15stage_locens_6 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_126:150")
two_strat_15stage_locens_7 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_51:75")
two_strat_15stage_locens_8 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_176:200")
two_strat_15stage_locens_9 <- read.csv("~/survrf/Outputs/15stage_500pt_locens_201:210")





pt500_15stage_low_observed <- c(two_strat_15stage_locens_1$observed, two_strat_15stage_locens_2$observed,
                                two_strat_15stage_locens_3$observed, two_strat_15stage_locens_4$observed,
                                two_strat_15stage_locens_5$observed, two_strat_15stage_locens_6$observed,
                                two_strat_15stage_locens_7$observed, two_strat_15stage_locens_8$observed,
                                two_strat_15stage_locens_9$observed)

pt500_15stage_low_IHsurvrf <- c(two_strat_15stage_locens_1$IHsurvrf, two_strat_15stage_locens_2$IHsurvrf,
                                two_strat_15stage_locens_3$IHsurvrf, two_strat_15stage_locens_4$IHsurvrf,
                                two_strat_15stage_locens_5$IHsurvrf, two_strat_15stage_locens_6$IHsurvrf,
                                two_strat_15stage_locens_7$IHsurvrf, two_strat_15stage_locens_8$IHsurvrf,
                                two_strat_15stage_locens_9$IHsurvrf)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_15stage_low_observed)), rep("IHsurvrf", length(pt500_15stage_low_IHsurvrf))),
  Value = c(pt500_15stage_low_observed, pt500_15stage_low_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_15stage_low_IHsurvrf)),
  design = rep("15 stages: Low Censoring Rate", 2*length(pt500_15stage_low_IHsurvrf)),
  setting = rep("2 Strata", 2*length(pt500_15stage_low_IHsurvrf))
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

lowcens_500pt_15stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(660, 780, by = 20), limits = c(660, 780)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_15stage_low_observed), mean(pt500_15stage_low_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_15stage_low_observed), mean(pt500_15stage_low_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



##################################################################
## plots for 500 pts, 15 stages, high censoring ##
## observational,                               ##
##################################################################

two_strat_15stage_hicens_1 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_101:125")
two_strat_15stage_hicens_2 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_126:150")
two_strat_15stage_hicens_3 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_1:25")
two_strat_15stage_hicens_4 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_26:50")
two_strat_15stage_hicens_5 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_51:75")
two_strat_15stage_hicens_6 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_151:175")
two_strat_15stage_hicens_7 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_76:100")
two_strat_15stage_hicens_8 <- read.csv("~/survrf/Outputs/15stage_500pt_highcens_176:200")
two_strat_15stage_hicens_9 <- read.csv("~/survrf/Outputs/15stage_500pt_hicens_201:205")


## using a threshold of 0.25
## 5 sims done
## 105 sims done
pt500_15stage_hi_observed <- c(two_strat_15stage_hicens_1$observed, two_strat_15stage_hicens_2$observed,
                               two_strat_15stage_hicens_3$observed, two_strat_15stage_hicens_4$observed,
                               two_strat_15stage_hicens_5$observed, two_strat_15stage_hicens_6$observed,
                               two_strat_15stage_hicens_7$observed, two_strat_15stage_hicens_8$observed,
                               two_strat_15stage_hicens_9$observed)
pt500_15stage_hi_IHsurvrf <- c(two_strat_15stage_hicens_1$IHsurvrf, two_strat_15stage_hicens_2$IHsurvrf,
                               two_strat_15stage_hicens_3$IHsurvrf, two_strat_15stage_hicens_4$IHsurvrf,
                               two_strat_15stage_hicens_5$IHsurvrf,two_strat_15stage_hicens_6$IHsurvrf,
                               two_strat_15stage_hicens_7$IHsurvrf, two_strat_15stage_hicens_8$IHsurvrf,
                               two_strat_15stage_hicens_9$IHsurvrf)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_15stage_hi_observed)), rep("IHsurvrf", length(pt500_15stage_hi_IHsurvrf))),
  Value = c(pt500_15stage_hi_observed, pt500_15stage_hi_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_15stage_hi_IHsurvrf)),
  design = rep("15 stages: High Censoring Rate", 2*length(pt500_15stage_hi_IHsurvrf)),
  setting = rep("2 Strata", 2*length(pt500_15stage_hi_IHsurvrf))
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

highcens_500pt_15stage_2strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(980, 1080, by = 20), limits = c(980, 1080)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_15stage_hi_observed), mean(pt500_15stage_hi_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_15stage_hi_observed), mean(pt500_15stage_hi_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 


##################################################################
## 1 strata: plots for 500 pts, 15 stages, low censoring ##
## observational,                               ##
##################################################################


one_strat_15stage_locens_1 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_1:25")
one_strat_15stage_locens_2 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_26:50")
one_strat_15stage_locens_3 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_101:125")
one_strat_15stage_locens_4 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_51:75")
one_strat_15stage_locens_5 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_76:100")
one_strat_15stage_locens_6 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_126:150")
one_strat_15stage_locens_7 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_151:175")
one_strat_15stage_locens_8 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_176:200")
one_strat_15stage_locens_9 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_locens_201:210")




pt500_15_lo1_observed <- c(one_strat_15stage_locens_1$observed, one_strat_15stage_locens_2$observed,
                           one_strat_15stage_locens_3$observed, one_strat_15stage_locens_4$observed,
                           one_strat_15stage_locens_5$observed, one_strat_15stage_locens_6$observed,
                           one_strat_15stage_locens_7$observed, one_strat_15stage_locens_8$observed,
                           one_strat_15stage_locens_9$observed)
pt500_15_lo1_IHsurvrf <- c(one_strat_15stage_locens_1$IHsurvrf, one_strat_15stage_locens_2$IHsurvrf,
                           one_strat_15stage_locens_3$IHsurvrf, one_strat_15stage_locens_4$IHsurvrf,
                           one_strat_15stage_locens_5$IHsurvrf, one_strat_15stage_locens_6$IHsurvrf,
                           one_strat_15stage_locens_7$IHsurvrf, one_strat_15stage_locens_8$IHsurvrf,
                           one_strat_15stage_locens_9$IHsurvrf)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_15_lo1_observed)), rep("IHsurvrf", length(pt500_15_lo1_IHsurvrf))),
  Value = c(pt500_15_lo1_observed, pt500_15_lo1_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_15_lo1_observed)),
  design = rep("15 stages: Low Censoring Rate", 2*length(pt500_15_lo1_observed)),
  setting = rep("1 Strata", 2*length(pt500_15_lo1_observed))
)


locens_500pt_15stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(660, 780, by = 20), limits = c(660, 780)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(. ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_15_lo1_observed), mean(pt500_15_lo1_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_15_lo1_observed), mean(pt500_15_lo1_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 
  



##################################################################
## 1 strata: plots for 500 pts, 15 stages, high censoring ##
## observational,                               ##
##################################################################

one_strat_15stage_hicens_1 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_51:75")
one_strat_15stage_hicens_2 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_1:25")
one_strat_15stage_hicens_3 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_26:50")
one_strat_15stage_hicens_4 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_76:100")
one_strat_15stage_hicens_5 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_151:175")
one_strat_15stage_hicens_6 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_101:125")
one_strat_15stage_hicens_7 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_126:150")
one_strat_15stage_hicens_8 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_highcens_176:200")
one_strat_15stage_hicens_9 <- read.csv("~/survrf/Outputs/1STRATA_15stage_500pt_hicens_201:210")


pt500_15_hi1_observed <- c(one_strat_15stage_hicens_1$observed, one_strat_15stage_hicens_2$observed,
                           one_strat_15stage_hicens_3$observed, one_strat_15stage_hicens_4$observed,
                           one_strat_15stage_hicens_5$observed, one_strat_15stage_hicens_6$observed,
                           one_strat_15stage_hicens_7$observed, one_strat_15stage_hicens_8$observed,
                           one_strat_15stage_hicens_9$observed)
pt500_15_hi1_IHsurvrf <- c(one_strat_15stage_hicens_1$IHsurvrf, one_strat_15stage_hicens_2$IHsurvrf,
                           one_strat_15stage_hicens_3$IHsurvrf, one_strat_15stage_hicens_4$IHsurvrf,
                           one_strat_15stage_hicens_5$IHsurvrf, one_strat_15stage_hicens_6$IHsurvrf,
                           one_strat_15stage_hicens_7$IHsurvrf, one_strat_15stage_hicens_8$IHsurvrf,
                           one_strat_15stage_hicens_9$IHsurvrf)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(pt500_15_hi1_observed)), rep("IHsurvrf", length(pt500_15_hi1_IHsurvrf))),
  Value = c(pt500_15_hi1_observed, pt500_15_hi1_IHsurvrf),
  n = rep("n = 500", 2*length(pt500_15_hi1_observed)),
  design = rep("15 stages: High Censoring Rate", 2*length(pt500_15_hi1_observed)),
  setting = rep("1 Strata", 2*length(pt500_15_hi1_observed))
)


hicens_500pt_15stage_1strata <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  scale_y_continuous(breaks = seq(980, 1080, by = 20), limits = c(980, 1080)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  facet_grid(setting ~ design + n) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(pt500_15_hi1_observed), mean(pt500_15_hi1_IHsurvrf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = -5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(pt500_15_hi1_observed), mean(pt500_15_hi1_IHsurvrf))),
             aes(x = Group, y = Value), color = "black", size = 3) 



####### grid plot

# Combine plots with a common label
plot_grid( lowcens_500pt_15stage_2strata + guides(fill = FALSE), 
          highcens_500pt_15stage_2strata + guides(fill = FALSE),
          locens_500pt_15stage_1strata,
          hicens_500pt_15stage_1strata, 
          align = "h", nrow = 2, ncol = 2)

# Save the plot with a higher resolution
ggsave("15stage.png", plot = last_plot(), dpi = 300, width = 10, height = 8)

## ---------------------------------------------------------------##
##        Equal censoring: 1 vs 2 strata -- lowcens               ##
## ---------------------------------------------------------------##


one_equal <- read.csv("~/survrf/Outputs/equal_1strata_10stage_locens")

one_observed <- c(one_equal$observed)

one_rf <- c(one_equal$IHsurvrf)


# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(one_observed)), rep("IHsurvrf", length(one_rf))),
  Value = c(one_observed, one_rf)
)

ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  #scale_y_continuous(breaks = seq(2200, 2500, by = 50), limits = c(2200, 2600)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(one_observed), mean(one_rf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = 5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(one_observed), mean(one_rf))),
             aes(x = Group, y = Value), color = "black", size = 3) + ggtitle("1 strata, 10 stages, 300 pts, low censoring, equal")


###### two strata equal

a <- read.delim2("~/survrf/Outputs/equal_2strata_10stage_locens.txt", header=FALSE, comment.char="#")
matching_indices_V1 <- which(apply(a, 1, function(row) any(grepl(".*no observed IHsurvrf*", row))))

selected_rows_V1 <- data.frame()

# Loop through each matching index
for (index in matching_indices_V1) {
  # Select the current row and the two rows after it, ensuring you don't exceed the data frame's bounds
  selected_rows_V1 <- rbind(selected_rows_V1, a[index:min(index+1, nrow(a)), ])
}

split_values_V1 <- strsplit(as.character(selected_rows_V1$X.1..4.479.7317.489.7977.483.5225.509.8735......NA......................0.), " ")
cleaned_split_values_V1 <- lapply(split_values_V1, function(x) x[x != ""])
sim_index <- sort(as.numeric(sapply(cleaned_split_values_V1, function(x) x[2])))

# Pull out numbers NOT in `sim_index`
not_in_sim_index <- setdiff(1:200, sim_index)

# Add 201 and include 149
result <- c(not_in_sim_index, 201, 149)

third_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[3])

fourth_values_V1 <- sapply(cleaned_split_values_V1, function(x) x[4])

two_rf <- as.numeric(fourth_values_V1)
two_observed <- as.numeric(third_values_V1)



# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(two_observed)), rep("IHsurvrf", length(two_rf))),
  Value = c(two_observed, two_rf)
)

ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  #scale_y_continuous(breaks = seq(2200, 2500, by = 50), limits = c(2200, 2600)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(two_observed), mean(two_rf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = 5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(two_observed), mean(two_rf))),
             aes(x = Group, y = Value), color = "black", size = 3) + ggtitle("2 strata, 10 stages, 300 pts, low censoring, equal")

####### 90% censoring, a = -6, and a = -3; 300 pts, 10 stages


a <- read.csv("~/survrf/Outputs/1STRATA_300pt_90_cens")

observed <- c(a$observed)

rf <- c(a$IHsurvrf)

sd(rf)

# Combine data into a dataframe
data <- data.frame(
  Group = c(rep("observed", length(observed)), rep("IHsurvrf", length(rf))),
  Value = c(observed, rf)
)

ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0, outlier.shape = NA, aes(color = Group)) +
  #scale_y_continuous(breaks = seq(2200, 2500, by = 50), limits = c(2200, 2600)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(color = Group), size = 3, alpha = 0.6) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank()) + theme(legend.position = "none") +
  geom_text(data = data.frame(Group = c("observed", "IHsurvrf"),
                              Value = c(mean(observed), mean(rf))),
            aes(x = Group, y = Value, label = round(Value, 2)), vjust = 5, hjust = 0,color = "black", size = 4) +
  geom_point(data = data.frame(Group = c("observed", "IHsurvrf"),
                               Value = c(mean(observed), mean(rf))),
             aes(x = Group, y = Value), color = "black", size = 3) + ggtitle("90% censoring, 10 stages, 300 pts")





