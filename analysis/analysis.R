# -----------------------------------------------------------------------------
# ANALYSIS OF PROPORTIONAL HAZARDS ASSUMPTION IN THE DIGIT-HF TRIAL
#
# This script performs a series of checks on the proportional hazards (PH)
# assumption for the Cox model applied to the DIGIT-HF trial data.
# The individual patient data (IPD) was reconstructed from the published
# Kaplan-Meier curves.
#
# Analysis steps:
# 0. Reproducing Cumulative Events Curves from DIGIT-HF trial
# 1. PH assumption check using scaled Schoenfeld residuals.
# 2. PH assumption check using log-log survival curves.
# 3. Modeling and plotting a time-dependent hazard ratio.
# 4. Modeling and plotting a Landmark Analysis at 12 months.
# -----------------------------------------------------------------------------


# 1. LOAD LIBRARIES AND DATA
# =============================================================================

# Use pacman to load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # For data manipulation and plotting
  broom,        # For tidy results
  survival,     # Core survival analysis functions
  survminer,    # For survival plots
  rstpm2,       # For flexible parametric survival models
  patchwork,    # To combine multiple plots
  broom,        # To tidy model outputs
  scales,       # For plot scales
  showtext      # For custom fonts
)

# Optional: Load custom font for aesthetics. 
# The script will run with system defaults if this fails.
tryCatch({
  font_add_google("Raleway", "raleway")
  showtext_auto()
}, warning = function(w) {
  message("Custom font 'Raleway' not found. Using system default.")
})


# Load the reconstructed Individual Patient Data (IPD)
# IMPORTANT: Assumes the CSV file is in the same directory as this script.
ipd_digit_hf <- read_csv("ipd_digit_hf1.csv")


# 2. FIT STANDARD COX MODEL AND TEST PH ASSUMPTION
# =============================================================================

# Fit the standard Cox Proportional Hazards model
digit_cox <- coxph(Surv(time, event) ~ arm, data = ipd_digit_hf)

# Verify the same HR and CI as the trial.
hr_ci <- tidy(digit_cox,exponentiate = T,conf.int = T) 

# --- 1. Create a formatted text string for the annotation ---
hr_text <- sprintf(
  "HR: %.2f (95%% CI: %.2f - %.2f)",
  hr_ci$estimate,
  hr_ci$conf.low,
  hr_ci$conf.high
)

# Perform the formal test for proportional hazards assumption
zph_test <- cox.zph(digit_cox)
zph_test

p_value <- zph_test$table[1, "p"]


# Plot 0:  Cumulative incidence of death from any cause or
#          first hospitalization for heart failure 
# -----------------------------------------------------------------------------
# Reproduction of trial curve for the main endpoint
cuminc_plot <- ggsurvplot(
  survfit(Surv(time, event) ~ arm, data = ipd_digit_hf),
  fun = "event",
  palette = c("#939598", "#026F99"),
  conf.int = T,
  censor = FALSE,
  legend.labs = c("Placebo", "Digitoxin"),
  legend.title = ""
)$plot +
  labs(
    title = "Death from Any Cause or First Hospitalization for Heart Failure",
    subtitle = "Replication of the DIGIT-HF trial curves",
    x = "Time in Months",
    y = "Cumulative Incidence"
  ) +
  annotate(
    "text", 
    x = 38,        # X-coordinate for the text
    y = 0.95,      # Y-coordinate for the text
    label = hr_text, 
    hjust = 0,     
    size = 4.5, 
    family = "raleway",
    fontface = "plain"
  ) +
  scale_x_continuous(
    limits = c(0, 110),                 
    breaks = scales::breaks_width(12)   
  ) +
  scale_y_continuous(
    limits = c(0, 1),                 
    breaks = scales::breaks_width(0.1), 
    labels = scales::percent_format()   
  ) +
  theme_minimal(base_family = "raleway") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, colour = "gray20"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 14, face = "plain"),
    legend.position = "top"
  )

cuminc_plot


# Plot 1: Scaled Schoenfeld Residuals vs. Time
# -----------------------------------------------------------------------------
# A non-zero slope in the smoothed line suggests a violation of the PH assumption.

# Extract residuals and corresponding times
residuals_df <- data.frame(
  residuals = residuals(digit_cox, type = "schoenfeld"),
  time = as.numeric(names(residuals(digit_cox, type = "schoenfeld")))
)

# Create the plot using ggplot
schoenfeld_plot <- ggplot(residuals_df, aes(x = time, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_point(alpha = 0.2, size = 2, color = "deepskyblue4") +
  # Using loess smoother to visualize the trend in residuals
  geom_smooth(method = "loess", formula = "y ~ x", 
              se = TRUE, color = "deepskyblue4", fill = "deepskyblue4", alpha = 0.15) +
  labs(
    title = "Schoenfeld Residuals vs. Time",
    subtitle = paste0("Formal Test for Non-Proportionality: p-value = ", format.pval(p_value, digits = 3, eps = 0.001)),
    x = "Time in Months",
    y = "Scaled Schoenfeld Residuals (arm)"
  ) +
  scale_x_continuous(breaks = breaks_width(12)) +
  theme_minimal(base_family = "raleway") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, colour = "gray20"),
    axis.title = element_text(size = 12)
  )

schoenfeld_plot


# Plot 2: Log-Log Survival Curves
# -----------------------------------------------------------------------------
# If the PH assumption holds, these curves should be parallel.
loglog_plot <- ggsurvplot(
  survfit(Surv(time, event) ~ arm, data = ipd_digit_hf),
  fun = "cloglog",
  palette = c("#939598", "#026F99"),
  conf.int = T,
  censor = FALSE,
  legend.labs = c("Placebo", "Digitoxin"),
  legend.title = ""
)$plot +
  labs(
    title = "Log-Log Survival Curves",
    subtitle = "Curves should be parallel if PH assumption holds",
    x = "log(Time in Months)",
    y = "log(-log(Survival Probability))"
  ) +
  theme_minimal(base_family = "raleway") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, colour = "gray20"),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

loglog_plot


# 3. MODEL AND PLOT THE TIME-DEPENDENT HAZARD RATIO
# =============================================================================

# Fit a flexible parametric model with a time-dependent effect for 'arm'
# This allows the Hazard Ratio to change over time.
fit_timecox <- stpm2(Surv(time, event) ~ arm, data = ipd_digit_hf, tvc = list(arm = 1))

# Generate predictions for the HR over time
times <- seq(min(ipd_digit_hf$time), max(ipd_digit_hf$time), length.out = 200)
pred_hr <- predict(
  fit_timecox,
  newdata = data.frame(arm = 1),
  grid = TRUE,
  full = TRUE,
  type = "hr",
  var = "arm",
  times = times,
  se.fit = TRUE
)


# Plot 3: Time-Dependent Hazard Ratio
# -----------------------------------------------------------------------------
# This plot shows how the treatment effect (HR) changes over the follow-up period.
time_dep_hr_plot <- ggplot(pred_hr, aes(x = time, y = Estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#026F99", alpha = 0.15) +
  geom_line(color = "#026F99", linewidth = 1.2) +
  geom_hline(yintercept = 1, color = "gray20", linetype = "dashed") +
  geom_hline(yintercept = exp(coef(digit_cox)), color = "deeppink4",
             linetype = "dotted", linewidth = 1) +
  annotate("text", x = 80, y = exp(coef(digit_cox))-0.1 + 0.04, 
           label = "Constant HR reported in Digit-HF trial (0.82)", 
           color = "deeppink4", family = "raleway") +
  scale_y_continuous(breaks = breaks_pretty(n = 8)) +
  scale_x_continuous(breaks = breaks_width(12)) +
  labs(
    title = "Time-Dependent Hazard Ratio of Digitoxin vs. Placebo",
    subtitle = "The effect of digitoxin appears to wane over time",
    x = "Time in Months",
    y = "Hazard Ratio (HR)"
  ) +
  theme_minimal(base_family = "raleway") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, colour = "gray20"),
    axis.title = element_text(size = 12)
  )

time_dep_hr_plot


# 4. COMBINE PLOTS AND SAVE
# =============================================================================

# Combine all three diagnostic plots into a single figure using patchwork
combined_plot <- (loglog_plot + schoenfeld_plot) / (time_dep_hr_plot) +
  plot_annotation(
    title = "Violation of Proportional Hazards Assumption in the DIGIT-HF Trial",
    caption = "Data reconstructed from published K-M curves",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", family = "raleway", hjust = 0.5),
      plot.caption = element_text(size = 10, family = "raleway", color = "gray40")
    )
  )

combined_plot


       # Save the final combined plot to a file
ggsave("ph_assumption_check_DIGIT-HF.png", 
       plot = combined_plot, 
       width = 12, 
       height = 10, 
       dpi = 300,
       bg = "white")


######## --- LANDMARK ANALYSIS: 12 MONTHS ---
# =============================================================================

# landmark= 12 months
landmark_time <- 12


# STEP 1: EARLY PERIOD ANALYSIS (0 to 12 months)
# -----------------------------------------------------------------------------

df_early <- ipd_digit_hf %>%
  mutate(
    time_early = if_else(time > landmark_time, landmark_time, time),
    event_early = if_else(time > landmark_time, 0, event) # 0 = censurado
  )

# Cox model for 12 months
cox_early <- coxph(Surv(time_early, event_early) ~ arm, data = df_early)

# Extract HR
hr_ci_early <- tidy(cox_early, exponentiate = TRUE, conf.int = TRUE)
hr_text_early <- sprintf(
  "HR for 0-12 months:\n%.2f (95%% CI: %.2f - %.2f)",
  hr_ci_early$estimate,
  hr_ci_early$conf.low,
  hr_ci_early$conf.high
)

# Early period plot
plot_early <- ggsurvplot(
  survfit(Surv(time_early, event_early) ~ arm, data = df_early),
  fun = "event",
  palette = c("#939598", "#026F99"),
  conf.int = FALSE,
  censor = FALSE,
  legend.labs = c("Placebo", "Digitoxin"),
  legend.title = "",
  xlim = c(0, landmark_time) 
)$plot +
  labs(
    title = "Early Period (0-12 Months)",
    x = "Time in Months",
    y = "Cumulative Incidence"
  ) +
  annotate("text", x = 1, y = 0.25, label = hr_text_early, hjust = 0, family = "raleway", fontface = "bold") +
  scale_y_continuous(limits = c(0, 0.3), labels = scales::percent) +
  theme_minimal(base_family = "raleway") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top", 
        legend.text = element_text(size=12))


# STEP 2: LATE PERIOD ANALYSIS (POST-12 months)
# -----------------------------------------------------------------------------

# Create late period dataset
df_landmark_cohort <- ipd_digit_hf %>%
  filter(time > landmark_time) %>%
  mutate(time_late = time - landmark_time)

# Model Cox for late period
cox_late <- coxph(Surv(time_late, event) ~ arm, data = df_landmark_cohort)

# Extract HR
hr_ci_late <- tidy(cox_late, exponentiate = TRUE, conf.int = TRUE)
hr_text_late <- sprintf(
  "HR for >12 months:\n%.2f (95%% CI: %.2f - %.2f)",
  hr_ci_late$estimate,
  hr_ci_late$conf.low,
  hr_ci_late$conf.high
)

# Late period plot
plot_late <- ggsurvplot(
  survfit(Surv(time_late, event) ~ arm, data = df_landmark_cohort),
  fun = "event",
  palette = c("#939598", "#026F99"),
  conf.int = FALSE,
  censor = FALSE,
  legend.title = "",
)$plot +
  labs(
    title = "Late Period (>12 Months)",
    x = "Time Since Landmark (Months)",
    y = ""
  ) +
  annotate("text", x = 5, y = 0.58, label = hr_text_late, hjust = 0, family = "raleway", fontface = "bold") +
  scale_y_continuous(limits = c(0, 0.7), labels = scales::percent) +
  theme_minimal(base_family = "raleway") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    #axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  )


# STEP 3: COMBINE PLOTS
# -----------------------------------------------------------------------------

combined_landmark_plot <- plot_early + plot_late +
  plot_annotation(
    title = "DIGIT-HF: Landmark Analysis at 12 Months for Primary Outcome",
    subtitle = "The treatment effect appears to be concentrated in the first year.",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "raleway"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, family = "raleway", color = "gray20")
    )
  )

combined_landmark_plot

# Optional: save the plot
# ggsave("landmark_analysis_12m.png", plot = combined_landmark_plot, width = 12, height = 7, dpi = 300, bg = "white")


