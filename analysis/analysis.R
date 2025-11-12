# -----------------------------------------------------------------------------
# ANALYSIS OF PROPORTIONAL HAZARDS ASSUMPTION IN THE DIGIT-HF TRIAL
#
# This script performs a series of checks on the proportional hazards (PH)
# assumption for the Cox model applied to the DIGIT-HF trial data.
# The individual patient data (IPD) was reconstructed from the published
# Kaplan-Meier curves.
#
# Analysis steps:
# 1. Visual check using scaled Schoenfeld residuals.
# 2. Visual check using log-log survival curves.
# 3. Modeling and plotting a time-dependent hazard ratio.
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

