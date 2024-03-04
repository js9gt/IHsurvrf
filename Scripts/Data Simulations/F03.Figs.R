
## code for plots

# Reshape data for ggplot2
## Plot 1, 2, 3
library(tidyr)

## Plot 1, 2, 3
library(dplyr)

# Plot 1, 2, 3
library(ggplot2)

# Writing the event.time summary as an excel file
library(writexl)




generate_plots_and_summary <- function(pts, num_patients, num_stages, other_text) {
  
  new_folder <- sprintf("~/Desktop/survrf/Figures/%spts_%sstage_%s", num_patients, num_stages, other_text)
  fs::dir_create(new_folder)  # Create a new folder
  
  # ------------------------------------------- #
  #          Plot 1: Event Time by Patient      #
  # ------------------------------------------- #
  
  event_times <- pts[, "event.time", ]
  max_stages <- ncol(event_times)
  
  p1 <- as.data.frame(event_times) %>%
    mutate(patient = rownames(event_times)) %>%
    pivot_longer(cols = 1:max_stages, names_to = "stage", values_to = "event.time", names_prefix = "V") %>%
    mutate(stage = as.numeric(stage),
           patient = as.numeric(gsub("pt", "", patient)))
  
  png(sprintf("%s/timeline_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p1_save <- ggplot(p1, aes(x = event.time, y = patient, color = factor(stage))) +
    geom_point() +
    geom_segment(aes(x = event.time, xend = event.time, y = patient, yend = patient + 0.25), size = 1) +
    geom_hline(yintercept = 1:max(p1$patient), linetype = "dashed", color = "gray") +
    scale_color_discrete() +
    labs(x = "Event Time", y = "Patient") +
    theme_minimal()
  
  print(p1_save)
  dev.off()
  
  # ------------------------------------------- #
  #              Plot 2: Average rates          #
  # ------------------------------------------- #
  
  rate_failure <- pts[, "rate.failure", ]
  rate_next_visit <- pts[, "rate.next.visit", ]
  rate_censoring <- pts[, "rate.censoring", ]
  
  p2 <- data.frame(
    stage = rep(1:max_stages, each = nrow(rate_failure)),
    rate_failure = c(rate_failure),
    rate_next_visit = c(rate_next_visit),
    rate_censoring = c(rate_censoring)
  )
  
  p2_avg <- p2 %>%
    group_by(stage) %>%
    summarize(
      avg_rate_failure = mean(rate_failure, na.rm = TRUE),
      avg_rate_next_visit = mean(rate_next_visit, na.rm = TRUE),
      avg_rate_censoring = mean(rate_censoring, na.rm = TRUE)
    )
  
  
  png(sprintf("%s/avg_rates_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p2_save <- ggplot(p2_avg, aes(x = stage, group = 1)) +
    geom_line(aes(y = avg_rate_failure, color = "Rate Failure"), size = 1) +
    geom_point(aes(y = avg_rate_failure, color = "Rate Failure"), size = 3) +
    geom_line(aes(y = avg_rate_next_visit, color = "Rate Next Visit"), size = 1) +
    geom_point(aes(y = avg_rate_next_visit, color = "Rate Next Visit"), size = 3) +
    geom_line(aes(y = avg_rate_censoring, color = "Rate Censoring"), size = 1) +
    geom_point(aes(y = avg_rate_censoring, color = "Rate Censoring"), size = 3) +
    labs(x = "Stage", y = "Average Value") +
    scale_color_manual(values = c("Rate Failure" = "blue", "Rate Next Visit" = "green", "Rate Censoring" = "red")) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(min(p2_avg$stage) - 1, max(p2_avg$stage), by = 1))
  
  print(p2_save)
  dev.off()
  
  # ------------------------------------------- #
  #             Plot 3: Pct Advancement         #
  # ------------------------------------------- #
  
  at_risk <- pts[, "at.risk", ]
  
  p3 <- data.frame(
    stage = rep(1:max_stages, each = nrow(at_risk)),
    at_risk = c(at_risk)
  )
  
  p3_percentage <- p3 %>%
    group_by(stage) %>%
    summarize(
      percentage_at_risk = sum(at_risk == 1) / n()
    )
  
  
  png(sprintf("%s/pct_atrisk_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p3_save <- ggplot(p3_percentage, aes(x = stage, y = percentage_at_risk, group = 1)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(x = "Stage", y = "Percentage of Patients at Risk") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(min(p3_percentage$stage) - 1, max(p3_percentage$stage), by = 5))
  
  print(p3_save)
  dev.off()
  
  # ------------------------------------------- #
  #              Plot 4: Average state          #
  # ------------------------------------------- #
  
  state_vals <- pts[, "state", ]
  
  p4 <- data.frame(
    stage = rep(1:max_stages, each = nrow(state_vals)),
    state_vals = c(state_vals)
  )
  
  p4_avg <- p4 %>%
    group_by(stage) %>%
    summarize(
      avg_state = mean(state_vals, na.rm = TRUE)
    )
  
  png(sprintf("%s/avg_states_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p4_save <- ggplot(p4_avg, aes(x = stage, group = 1)) +
    geom_line(aes(y = avg_state, color = "Average State Value"), size = 1) +
    geom_point(aes(y = avg_state, color = "Average State Value"), size = 3) +
    labs(x = "Stage", y = "Average Value") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(min(p4_avg$stage)-1, max(p4_avg$stage), by = 5))
  
  print(p4_save)
  dev.off()
  
  # ------------------------------------------- #
  #               Plot 5: Pct Action            #
  # ------------------------------------------- #
  
  action_vals <- pts[, "action", ]
  
  p5 <- data.frame(
    stage = rep(1:max_stages, each = nrow(action_vals)),
    action_vals = c(action_vals)
  )
  
  p5_percentage <- p5 %>%
    group_by(stage) %>% 
    summarize(percentage = mean(action_vals == 1, na.rm = TRUE) * 100,
              denominator_count = sum(!is.na(action_vals)))
  

  png(sprintf("%s/pct_action_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p5_save <- ggplot(p5_percentage, aes(x = stage, y = percentage)) +
    geom_line(aes(y = percentage), size = 1) +
    geom_point(aes(y = percentage), size = 3) +
    geom_text(aes(label = paste("n =", denominator_count)), 
              vjust = -1, hjust = 0.5, size = 3, color = "darkred") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "blue") +
    theme_minimal() + 
    labs(x = "Stage", y = "Percentage of Action = 1") +
    scale_x_continuous(breaks = seq(min(p5_percentage$stage) - 1, max(p5_percentage$stage), by = 5))
  
  print(p5_save)
  dev.off()
  
  # ------------------------------------------- #
  #            Plot 6: Pct Gamma Delta          #
  # ------------------------------------------- #
  
  gamma <- pts[, "gamma", ]
  delta <- pts[, "delta", ]
  
  p6 <- data.frame(
    stage = rep(1:max_stages, each = nrow(gamma)),
    gamma = c(gamma),
    delta = c(delta)
  )
  
  p6_percentage <- p6 %>%
    group_by(stage) %>% 
    summarize(pct_gamma = mean(gamma == 1, na.rm = TRUE) * 100,
              pct_censor = mean(delta == 0, na.rm = TRUE) * 100,
              gamma_count = sum(!is.na(gamma)),
              censor_count = sum(!is.na(delta)))
  
  
  png(sprintf("%s/pct_gamma_delta_%spts_%sstage_%s.png", new_folder, num_patients, num_stages, other_text), units = "px", res = 75)
  
  p6_save <- ggplot(p6_percentage, aes(x = stage, group = 1)) +
    geom_line(aes(y = pct_gamma, color = "Percentage of Failures"), size = 1) +
    geom_point(aes(y = pct_gamma, color = "Percentage of Failures"), size = 3) +
    geom_line(aes(y = pct_censor, color = "Percentage of Censoring"), size = 1) +
    geom_point(aes(y = pct_censor, color = "Percentage of Censoring"), size = 3) +
    labs(x = "Stage", y = "Percentage") + theme_minimal()  +
    scale_x_continuous(breaks = seq(min(p6_percentage$stage) - 1, max(p6_percentage$stage), by = 5))
  
  print(p6_save)
  dev.off()
  
  # ------------------------------------------- #
  #            Stage Event Time Summary         #
  # ------------------------------------------- #
  
  event_time_vals <- pts[, "event.time", , drop = FALSE]
  event_time_long <- as.data.frame(as.table(event_time_vals))
  
  summary_table <- event_time_long %>%
    group_by(Var3) %>%
    summarise(
      Mean = mean(Freq, na.rm = TRUE),
      Median = median(Freq, na.rm = TRUE),
      Min = min(Freq, na.rm = TRUE),
      Max = max(Freq, na.rm = TRUE),
      N = sum(!is.na(Freq))
    ) %>%
    rename(Stage = Var3)
  
  write_xlsx(summary_table, path =sprintf("~/Desktop/survrf/Outputs/event.time_summary_%spts_%sstage_%s.xlsx", num_patients, num_stages, other_text))
  
  
  message("Plots and summary table generated successfully!")
}


# Example usage:
# generate_plots_and_summary(pts)



