# ========================================================================================================
# R Script: Post-Anthesis Rice Grain Filling & iWUE Dynamics
# Author: Venice Margarette B. Juanillas
# Vibe coding done with Google Gemini
# Location: IRRI Zeigler Experiment Station, Los Baños, Philippines (14.167 N, 121.242 E)
# Description: This R script features: 6-ODE System, Multi-Scenario Sweeps, Multi-Output OAT & Morris GSA
#                   (Evaluates 13 Parameters + 2 Initial States across 5 Outputs)
# ========================================================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
required_packages <- c("deSolve", "tidyverse", "patchwork", "nasapower", "ggrepel", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(deSolve)
library(tidyverse)
library(patchwork)
library(nasapower)
library(ggrepel)

# ------------------------------------------------------------------------------
# 2. Fetch & Parse Microclimate Data (NASA POWER)
# ------------------------------------------------------------------------------
cat("Fetching NASA POWER weather data for IRRI Zeigler Station...\n")

# can we get other environmental drivers? VPD perhaps?
irri_weather_raw <- get_power(
  community = "ag",
  pars = c("T2M_MAX", "T2M_MIN", "ALLSKY_SFC_SW_DWN"),
  lonlat = c(121.242, 14.167),
  dates = c("2025-07-16", "2025-08-16"),
  temporal_api = "daily"
)

irri_df <- as.data.frame(irri_weather_raw)
date_col <- intersect(c("YYYYMMDD", "YYYYMMDD_DATE", "dates", "DATE"), names(irri_df))[1]

irri_weather <- irri_df %>%
  mutate(date_clean = as.Date(as.character(.data[[date_col]]), format = "%Y%m%d")) %>%
  filter(!is.na(ALLSKY_SFC_SW_DWN), ALLSKY_SFC_SW_DWN > -900) %>%
  mutate(
    day_idx  = row_number() - 1,
    PAR_peak = (ALLSKY_SFC_SW_DWN * 0.45 * 4.56 * 1e6 / (12 * 3600)) * (pi / 2.0)
  )

tmin_interp <- approxfun(irri_weather$day_idx, irri_weather$T2M_MIN, rule = 2)
tmax_interp <- approxfun(irri_weather$day_idx, irri_weather$T2M_MAX, rule = 2)
par_interp  <- approxfun(irri_weather$day_idx, irri_weather$PAR_peak, rule = 2)

get_irri_environment <- function(t, T_night_add = 0, T_max_add = 0, I_scale = 1.0) {
  day    <- floor(t)
  hour   <- (t - day) * 24.0
  T_min  <- tmin_interp(t) + T_night_add
  T_max  <- tmax_interp(t) + T_max_add
  I_peak <- par_interp(t) * I_scale
  
  T_mean <- (T_max + T_min) / 2.0
  T_amp  <- (T_max - T_min) / 2.0
  T_curr <- T_mean + T_amp * sin((hour - 8.0) * (2.0 * pi / 24.0))
  
  I_curr <- if (hour >= 6.0 && hour <= 18.0) max(0.0, I_peak * sin((hour - 6.0) * (pi / 12.0))) else 0.0
  return(list(T_curr = T_curr, I_curr = I_curr, T_max_day = T_max))
}

# ------------------------------------------------------------------------------
# 3. Model Parameters & 6-ODE Dynamic System
# ------------------------------------------------------------------------------
base_params <- c(
  g_s_max     = 0.60,   K_I         = 300.0,  A_max_g     = 25.0,
  g_0         = 0.08,   R_leaf_0    = 1.50,   Q10         = 2.0,
  M           = 30.0,   eta_starch  = 0.83,   k_remob     = 0.08,
  D_grain     = 18.0,   mu_sen      = 0.035,  beta_stress = 1.50,
  g_c         = 0.005,  dry_rate    = 0.020,  T_night_add = 0.0,
  T_max_add   = 0.0,    I_scale     = 1.0
)

base_y0 <- c(W_grain = 0.0, NSC_stem = 200.0, LAI = 5.0)
t_end   <- max(irri_weather$day_idx)
times   <- seq(0, t_end, by = 0.01)

rice_irri_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    env       <- get_irri_environment(t, T_night_add, T_max_add, I_scale)
    T_curr    <- env$T_curr
    I_curr    <- env$I_curr
    T_max_day <- env$T_max_day
    
    # 1. Soil Water & Stomatal Conductance (Eq 1)
    S_w    <- max(0.05, 1.0 - dry_rate * t)
    f_Tair <- exp(-((T_curr - 28.0) / 8.0)^2)
    g_s    <- g_s_max * (I_curr / (I_curr + K_I)) * f_Tair * S_w
    
    # 2. Respiration & Net Photosynthesis (Eq 2)
    R_leaf <- R_leaf_0 * (Q10^((T_curr - 20.0) / 10.0))
    if (I_curr > 5.0) {
      A_gross <- A_max_g * (g_s / (g_s + g_0))
      A_net   <- max(0.0, A_gross - R_leaf)
    } else {
      A_net   <- -R_leaf
    }
    
    # 3. Intrinsic Water-Use Efficiency (Eq 3)
    WUE_i <- if (I_curr > 5.0 && g_s > 0.05) A_net / (g_s + g_c) else 0.0
    
    # 4. Canopy NPP & Direct Supply (Eq 4)
    NPP_instant   <- 2.592 * LAI * A_net
    direct_supply <- eta_starch * max(0.0, NPP_instant)
    
    # 5. Indirect Stem Reserve Remobilization (Eq 5)
    deficit <- max(0.0, D_grain - direct_supply)
    F_remob <- min(k_remob * max(0.0, NSC_stem), deficit)
    dNSC_dt <- -F_remob
    
    # 6. Grain Biomass Accumulation & Leaf Senescence (Eq 6)
    f_heat  <- max(0.0, 1.0 - max(0.0, T_max_day - 35.0) * 0.05)
    dW_dt   <- min(D_grain, direct_supply + F_remob) * f_heat
    dLAI_dt <- -mu_sen * LAI * (1.0 + beta_stress * (1.0 - S_w))
    
    return(list(
      c(dW_dt, dNSC_dt, dLAI_dt),
      g_s = g_s, A_net = A_net, WUE_i = WUE_i, NPP_instant = NPP_instant
    ))
  })
}

# Helper to run simulation and summarize all 5 state/flux outputs
run_sim_summary <- function(custom_params = list(), custom_y0 = base_y0) {
  p  <- base_params
  y0 <- custom_y0
  
  # Route initial state overrides (NSC_stem_0, LAI_0) vs kinetic parameters
  for (nm in names(custom_params)) {
    if (nm == "NSC_stem_0") {
      y0["NSC_stem"] <- custom_params[[nm]]
    } else if (nm == "LAI_0") {
      y0["LAI"] <- custom_params[[nm]]
    } else {
      p[nm] <- custom_params[[nm]]
    }
  }
  
  out <- as.data.frame(ode(y = y0, times = times, func = rice_irri_ode, parms = p, method = "lsoda"))
  
  daytime_df <- out %>%
    mutate(hour_of_day = (time - floor(time)) * 24.0) %>%
    filter(hour_of_day >= 7.0 & hour_of_day <= 17.0 & A_net > 0 & g_s > 0.05)
  
  dt <- diff(out$time)[1]
  
  c(
    W_grain   = tail(out$W_grain, 1),
    Cum_NPP   = sum(out$NPP_instant) * dt,
    Mean_WUEi = mean(daytime_df$WUE_i, na.rm = TRUE),
    Mean_LAI  = mean(out$LAI, na.rm = TRUE),
    Final_NSC = tail(out$NSC_stem, 1)
  )
}

# ------------------------------------------------------------------------------
# 4. Stress Scenario Configurations
# ------------------------------------------------------------------------------
scenarios <- list(
  "Baseline"        = list(params = list(), y0 = base_y0),
  "Low Stem NSC"    = list(params = list(), y0 = c(W_grain = 0, NSC_stem = 100, LAI = 5)),
  "Extreme Drought" = list(params = list(dry_rate = 0.08, beta_stress = 2.5), y0 = base_y0),
  "Extreme Heat"    = list(params = list(T_max_add = 5.0, T_night_add = 6.0), y0 = base_y0),
  "Monsoon Shade"   = list(params = list(I_scale = 0.35), y0 = base_y0)
)

# Combined list of 13 physiological parameters + 2 initial states (NSC_stem_0, LAI_0)
all_inputs <- c(
  "D_grain", "k_remob", "NSC_stem_0", "LAI_0", "A_max_g",
  "g_s_max", "g_0", "K_I", "R_leaf_0", "Q10",
  "eta_starch", "mu_sen", "beta_stress", "dry_rate", "g_c"
)
output_vars <- c("W_grain", "Cum_NPP", "Mean_WUEi", "Mean_LAI", "Final_NSC")

# Helper to retrieve scenario-specific baseline value for any parameter or initial state
get_scen_base_val <- function(inp_name, scen_cfg) {
  if (!is.null(scen_cfg$params[[inp_name]])) return(scen_cfg$params[[inp_name]])
  if (inp_name == "NSC_stem_0") return(scen_cfg$y0[["NSC_stem"]])
  if (inp_name == "LAI_0")      return(scen_cfg$y0[["LAI"]])
  return(base_params[[inp_name]])
}

# ------------------------------------------------------------------------------
# 5. One-at-a-Time (OAT) Sensitivity Analysis Module (Parameters + Initial States)
# ------------------------------------------------------------------------------
cat("Executing One-at-a-Time (OAT) Sensitivity Analysis across all scenarios...\n")

oat_results <- list()

for (scen_name in names(scenarios)) {
  scen_cfg <- scenarios[[scen_name]]
  base_out <- run_sim_summary(scen_cfg$params, scen_cfg$y0)
  
  for (p_name in all_inputs) {
    base_val <- get_scen_base_val(p_name, scen_cfg)
    
    # -20% Evaluation
    p_low <- scen_cfg$params
    p_low[[p_name]] <- base_val * 0.80
    out_low <- run_sim_summary(p_low, scen_cfg$y0)
    
    # +20% Evaluation
    p_high <- scen_cfg$params
    p_high[[p_name]] <- base_val * 1.20
    out_high <- run_sim_summary(p_high, scen_cfg$y0)
    
    for (ov in output_vars) {
      b_y <- base_out[[ov]]
      pct_low  <- if (abs(b_y) > 1e-6) ((out_low[[ov]] - b_y) / b_y) * 100 else 0
      pct_high <- if (abs(b_y) > 1e-6) ((out_high[[ov]] - b_y) / b_y) * 100 else 0
      elasticity <- ((pct_low / -20) + (pct_high / 20)) / 2.0
      
      oat_results[[length(oat_results) + 1]] <- data.frame(
        Scenario   = scen_name,
        Input      = p_name,
        Output     = ov,
        Base_Val   = base_val,
        Out_Low    = out_low[[ov]],
        Out_High   = out_high[[ov]],
        Pct_Low    = pct_low,
        Pct_High   = pct_high,
        Elasticity = elasticity
      )
    }
  }
}

oat_df <- bind_rows(oat_results)
write.csv(oat_df, "OAT_Sensitivity_All_Scenarios_With_InitStates.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Global Sensitivity Analysis (Morris Screening) Module
# ------------------------------------------------------------------------------
cat("Executing Morris Global Sensitivity Analysis (GSA) including initial states...\n")

r_trajectories <- 15
p_levels <- 4
delta <- p_levels / (2 * (p_levels - 1))
k_vars <- length(all_inputs)

gsa_results <- list()

set.seed(42)
for (scen_name in names(scenarios)) {
  scen_cfg <- scenarios[[scen_name]]
  
  # Bounds (+/- 25%) anchored to each scenario's true baseline
  b_min <- sapply(all_inputs, function(p) get_scen_base_val(p, scen_cfg) * 0.75)
  b_max <- sapply(all_inputs, function(p) get_scen_base_val(p, scen_cfg) * 1.25)
  
  ee_list <- lapply(output_vars, function(x) matrix(0, nrow = r_trajectories, ncol = k_vars))
  names(ee_list) <- output_vars
  
  for (i in 1:r_trajectories) {
    x_base_grid <- floor(runif(k_vars) * p_levels) / (p_levels - 1)
    p_phys_base <- b_min + x_base_grid * (b_max - b_min)
    
    p_base_full <- scen_cfg$params
    p_base_full[all_inputs] <- p_phys_base
    y_base <- run_sim_summary(p_base_full, scen_cfg$y0)
    
    for (j in 1:k_vars) {
      x_pert <- x_base_grid
      x_pert[j] <- if (x_base_grid[j] + delta <= 1) x_base_grid[j] + delta else x_base_grid[j] - delta
      
      p_phys_pert <- b_min + x_pert * (b_max - b_min)
      p_pert_full <- scen_cfg$params
      p_pert_full[all_inputs] <- p_phys_pert
      y_pert <- run_sim_summary(p_pert_full, scen_cfg$y0)
      
      delta_actual <- (p_phys_pert[j] - p_phys_base[j]) / (b_max[j] - b_min[j])
      for (ov in output_vars) {
        ee_list[[ov]][i, j] <- if (delta_actual != 0) (y_pert[[ov]] - y_base[[ov]]) / delta_actual else 0
      }
    }
  }
  
  for (ov in output_vars) {
    gsa_results[[length(gsa_results) + 1]] <- data.frame(
      Scenario = scen_name,
      Input    = all_inputs,
      Output   = ov,
      Mu_Star  = colMeans(abs(ee_list[[ov]])),
      Sigma    = apply(ee_list[[ov]], 2, sd)
    )
  }
}

gsa_df <- bind_rows(gsa_results) %>%
  group_by(Scenario, Output) %>%
  mutate(Norm_Mu_Star = if (max(Mu_Star) > 0) Mu_Star / max(Mu_Star) else 0) %>%
  ungroup()

write.csv(gsa_df, "GSA_Morris_All_Scenarios_With_InitStates.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Visualization: Tornado Plots, Multi-Output Heatmaps & GSA Scatter Plots
# ------------------------------------------------------------------------------


# Daily Environmental Drivers Summary (Baseline Scenario)
daily_env_drivers <- irri_weather %>%
  group_by(calendar_date) %>%
  summarise(
    PAR_peak = max(I_curr, na.rm = TRUE),
    T_max    = max(T_curr, na.rm = TRUE),
    T_mean   = mean(T_curr, na.rm = TRUE),
    T_min    = min(T_curr, na.rm = TRUE),
    .groups = "drop"
  )

# Panel 1: Driving Solar Irradiance (PAR Peak)
p_par <- ggplot(daily_env_drivers, aes(x = calendar_date, y = PAR_peak)) +
  geom_area(fill = "#FFF3E0", alpha = 0.6) +
  geom_line(color = "#E65100", linewidth = 1.1) +
  geom_point(color = "#BF360C", size = 1.8) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  labs(title = "PAR Irradiance Peak (I_peak)",
       x = "Date", y = expression(I[peak] ~ (mu * mol ~ m^-2 ~ s^-1)))

# Panel 2: Canopy Temperature Dynamics (T_max, T_mean, T_min)
p_temp <- ggplot(daily_env_drivers, aes(x = calendar_date)) +
  geom_ribbon(aes(ymin = T_min, ymax = T_max), fill = "#E1F5FE", alpha = 0.5) +
  geom_line(aes(y = T_max, color = "T_max"), linewidth = 0.9) +
  geom_line(aes(y = T_mean, color = "T_mean"), linewidth = 1.1) +
  geom_line(aes(y = T_min, color = "T_min"), linewidth = 0.9) +
  scale_color_manual(
    name = NULL,
    values = c("T_max" = "#D32F2F", "T_mean" = "#1976D2", "T_min" = "#388E3C"),
    labels = c("T_max" = "Daily Max", "T_mean" = "Daily Mean", "T_min" = "Daily Min")
  ) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 days") +
  labs(title = "Diurnal Canopy Temperature Dynamics (T_curr)",
       x = "Date", y = expression("Temperature (" * degree * "C)")) +
  theme(legend.position = "top")

# Plot environmetal Factors first
env_panel <- (p_par | p_temp) +
  plot_annotation(
    title    = "Environmental Factors",
    subtitle = "Solar Irradiance and Temperature at IRRI Station",
    caption  = "Source: NASA  Power (July 16 – August 16, 2025)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 10)),
      plot.caption  = element_text(size = 9, face = "italic", color = "gray40")
    )
  ) &
  plot_layout(guides = "collect") &
  theme(legend.position = "top")



# Helper function to generate horizontal Tornado plots for any model output
make_tornado_plot <- function(df_oat, target_output, plot_title, x_axis_label) {
  df_oat %>%
    filter(Output == target_output) %>%
    group_by(Scenario) %>%
    mutate(Input_Ordered = reorder(Input, abs(Elasticity))) %>%
    ungroup() %>%
    ggplot(aes(y = Input_Ordered)) +
    geom_bar(aes(x = Pct_High, fill = "+20% Perturbation"), stat = "identity", width = 0.65, color = "black", linewidth = 0.2) +
    geom_bar(aes(x = Pct_Low, fill = "-20% Perturbation"), stat = "identity", width = 0.65, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
    facet_wrap(~Scenario, scales = "free") +
    scale_fill_manual(values = c("-20% Perturbation" = "#D84315", "+20% Perturbation" = "#1565C0")) +
    labs(
      title    = plot_title,
      subtitle = "Percentage change under ±20% parameter and initial state perturbations",
      x        = x_axis_label,
      y        = "Model Parameter / Initial State",
      fill     = "Perturbation"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 10.5),
      plot.title       = element_text(face = "bold", size = 13)
    )
}

# 1. Tornado Plot: Final Grain Weight (W_grain)
p_tornado_wgrain <- make_tornado_plot(
  oat_df,
  target_output = "W_grain",
  plot_title    = "Final Grain Weight (W_grain) Tornado Sensitivity Plot",
  x_axis_label  = "Percentage Change in Final Grain Weight (%)"
)

# 2. Tornado Plot: Cumulative Canopy NPP (Cum_NPP)
p_tornado_npp <- make_tornado_plot(
  oat_df,
  target_output = "Cum_NPP",
  plot_title    = "Cumulative Canopy NPP Tornado Sensitivity Plot",
  x_axis_label  = "Percentage Change in Cumulative Canopy NPP (%)"
)

# 3. Tornado Plot: Mean Daytime Intrinsic WUE (Mean_WUEi)
p_tornado_wuei <- make_tornado_plot(
  oat_df,
  target_output = "Mean_WUEi",
  plot_title    = "Intrinsic Water-Use Efficiency (WUE_i) Tornado Sensitivity Plot",
  x_axis_label  = "Percentage Change in Mean Daytime WUE_i (%)"
)

# 4. Multi-Output OAT Elasticity Heatmap Across All Scenarios
p_oat_heatmap <- ggplot(oat_df, aes(x = Output, y = factor(Input, levels = rev(all_inputs)), fill = Elasticity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Elasticity)), size = 2.5) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828", midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish) +
  facet_wrap(~Scenario, nrow = 1) +
  labs(
    title = "OAT Sensitivity Elasticity (13 Parameters + 2 Initial States vs. 5 Model Outputs)",
    x = "Model Output Variable", y = "Parameter / Initial State", fill = "Elasticity"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"), axis.text.y = element_text(face = "bold"))

# 5. Morris Global Sensitivity Heatmap Across All Scenarios
p_gsa_heatmap <- ggplot(gsa_df, aes(x = Output, y = factor(Input, levels = rev(all_inputs)), fill = Norm_Mu_Star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Norm_Mu_Star)), size = 2.5) +
  scale_fill_gradient(low = "#FFFDE7", high = "#0D47A1", limits = c(0, 1)) +
  facet_wrap(~Scenario, nrow = 1) +
  labs(
    title = "Morris Global Sensitivity (Normalized μ* Across All Parameters, Initial States & Outputs)",
    x = "Model Output Variable", y = "Parameter / Initial State", fill = "Rel. μ*"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"), axis.text.y = element_text(face = "bold"))

# 6. Morris Screening Scatter Plot (Mu* vs Sigma) on Final Grain Weight (D_grain unclipped)
p_gsa_yield <- gsa_df %>%
  filter(Output == "W_grain") %>%
  ggplot(aes(x = Mu_Star, y = Sigma, color = Input, label = Input)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  facet_wrap(~Scenario, scales = "free") +
  labs(
    title = "Morris Screening (μ* vs σ) on Final Grain Weight (W_grain) — Including Initial States",
    x = expression(bold("Mean Absolute Effect (" * mu * "* - Overall Influence)")),
    y = expression(bold("Standard Deviation (" * sigma * " - Non-linearity / Interactions)"))
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Export all Tornado plots, Heatmaps, and Morris GSA Scatter plots into a single multi-page PDF
pdf("Rice_Model_Full_Sensitivity_With_Tornado.pdf", width = 15, height = 9)
print(p_tornado_wgrain)
print(p_tornado_npp)
print(p_tornado_wuei)
print(p_oat_heatmap)
print(p_gsa_heatmap)
print(p_gsa_yield)
dev.off()

cat("Master R Script finished successfully. Generated outputs:
- OAT_Sensitivity_All_Scenarios_With_InitStates.csv
- GSA_Morris_All_Scenarios_With_InitStates.csv
- Rice_Model_Full_Sensitivity_With_Tornado.pdf\n")
