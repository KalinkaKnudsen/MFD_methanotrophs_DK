# Script for importing and displaying GC-MS for the 13CH4 incubations
# Kalinka Sand Knudsen

setwd("path/to/your/repo/MFD_methanotrophs_DK/")


library(tidyverse)
library(vroom)
library(ggplot2)
library(readxl)
library(lubridate)

# Path to your Excel file
excel_file <- "data/oxi_rate_raw_data.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(excel_file)

# Read and combine all sheets
df_all <- map_df(sheet_names, ~ read_excel(excel_file, sheet = .x) %>%
                   mutate(sheet_name = .x))  # optional: track origin

# Filter out failed injections, standardize a calibration Type, parse injection timestamps,
# and compute `oxi_time` (minutes since first injection per `sample_id`).
df_all <- df_all %>%
  filter(!`Injection Name`=="7_fail")%>%
  mutate(Type=if_else(Type=="Unknown" & `Injection Name`=="2_5 ppm CH4 CO H2", "Calibration_check", Type))%>%
  mutate(
    Inject_Time = ymd_hms(`Inject Time`, tz = "UTC"),
    sample_id = if_else(Type=="Unknown",sub("_.*", "", `Injection Name`), NA)
  ) %>%
  group_by(sample_id) %>%
  mutate(
    oxi_time = as.numeric(difftime(Inject_Time, min(Inject_Time), units = "mins"))
  ) %>%
  ungroup()


# Read sample metadata / linkage table and join to `df_all` to attach dates, masses, etc.
sample_link<-readxl::read_excel("data/Sample_linkage_mass.xlsx")%>%
  mutate(sample_id=as.character(sample_id))%>%mutate(Date=as.Date(Date))%>%filter(!Date=="2025-04-22")%>%distinct()

# Left-join sample metadata and set `oxi_time` to NA for calibration standards
df_all<-df_all%>%left_join(sample_link)%>%mutate(oxi_time=if_else(Type=="Calibration Standard", NA, oxi_time))



# Quick check plot: CH4 concentration over time per `sample_id`
df_all%>%filter(Component=="CH4 [FID]")%>%
  ggplot(., aes(x = oxi_time, y = Amount, colour=sample_id)) +
  geom_line() + # Use geom_line() for a line plot
  labs(title = "CH4 Concentration Over Time",
    x = "Date and Time",
    y = "CH4 [ppm]") +
  theme_minimal()




############### Okay the data is now looking fine ###################
## Next up,  1) plot the data as point for Methylocella only and 2) calculate the rate using 1. order kinetics 

my_pub_theme <- theme_classic() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    # Axis titles & text
    axis.title.x = element_text(size = 5, margin = margin(t = 0.1, b=1)),
    axis.title.y = element_text(size = 5, margin = margin(r = 0.05)),
    axis.text = element_text(size = 5),
    
    # Force axis lines to be drawn
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    
    # Plot title & subtitle
    plot.subtitle = element_text(
      size = 5,
      #  hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    plot.title = element_text(
      size = 5,
      face="bold",
      #  hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    
    # Legend text & title
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    
    # Legend spacing & key size
    legend.key.size = unit(0.1, "cm"),       # overall key box size
    # legend.key.height = unit(0.2, "cm"),     # vertical spacing
    #  legend.key.width = unit(0.2, "cm"),      # horizontal spacing
    legend.spacing.y = unit(0.07, "cm"),      # space between rows
    legend.spacing.x = unit(0.07, "cm"),      # space between columns
    legend.position = c(0.8, 0.8),
    legend.justification = c(0, 1),  
    # Legend background
    legend.background = element_rect(fill = "transparent")
  )


custom_colors <- c(
  "10" =   "#E41A1C", # f7081d
  "11" =  "brown", # Purple
  "12" =  "#FC8D62", # Coral
  "9" =   "orange",   # gold
  
  "5" ="#7677cd", # Teal
  "6" = "blue4", # Periwinkle
  "7" ="#2219e6", # Blue
  "8" ="#4db2a0", 
  "1" = "#006400",  # DarkGreen
  "2" = "#228B22",  # ForestGreen
  "3" = "#32CD32",  # LimeGreen
  "4" = "#7CFC00"   # LawnGreen
  
)


df_all%>%filter(Component=="CH4 [FID]")%>%
  filter(location=="Fussingø")%>%
  ggplot(aes(x = oxi_time, y = Amount)) +
  geom_point(aes(fill = sample_id), size = 1, alpha = 0.6,
             color = "black", pch = 21, stroke = 0.3) +
  scale_fill_manual(values = custom_colors)+ #, name = paste0(unique(df_all$sample_full[df_all$location == "Fussingø"]), " – Fussingø"))+
  labs(x = "Time [min]",
       y = expression(CH[4]~"[ppm]"),
       subtitle = paste0(unique(df_all$sample_full[df_all$location == "Fussingø"]), " – Fussingø"))+
  my_pub_theme +
  theme(legend.position = "none")


# Filter data for one sample (e.g., sample_id == "10")
df_sample <- df_all %>%
  filter(Component == "CH4 [FID]", sample_id == "10")

lm_fit <- lm(log(Amount) ~ oxi_time, data = df_sample)


# Extract slope (negative rate constant)
k <- -coef(lm_fit)[2]
C0 <- df_sample$Amount[1]
df_sample <- df_sample %>%
  mutate(predicted = C0 * exp(-k * oxi_time))%>%arrange(oxi_time)

time_seq <- seq(min(df_sample$oxi_time), max(df_sample$oxi_time), length.out = 200)
decay_curve <- data.frame(
  oxi_time = time_seq,
  CH4_pred = C0 * exp(-k * time_seq)
)

ggplot(df_sample, aes(x = oxi_time, y = Amount)) +
  geom_point(size = 1.5, fill = "skyblue", color = "black", pch = 21) +
  geom_line(data = decay_curve, aes(x = oxi_time, y = CH4_pred),
            color = "darkred", linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Time [min]",
    y = expression(CH[4]~"[ppm]"),
    title = paste("First-order decay (k =", round(k, 4), ")")
  ) +
  theme_classic()


# Filter CH4 data
df_ch4 <- df_all %>%
  filter(Component == "CH4 [FID]") %>%
  filter(!is.na(location))%>%
  filter(oxi_time < 600) %>%
  arrange(sample_id, oxi_time)

# Initialize empty lists
curve_list <- list()
data_list <- list()
meta_list <- list()

# Loop over each sample
for (id in unique(df_ch4$sample_id)) {
  df_sample <- df_ch4 %>% filter(sample_id == id)
  
  # Skip if not enough points
  if (nrow(df_sample) < 3) next
  
  # Fit log-linear model
  lm_fit <- lm(log(Amount) ~ oxi_time, data = df_sample)
  k <- -coef(lm_fit)[2]
  C0 <- df_sample$Amount[1]
  
  # Generate decay curve
  time_seq <- seq(min(df_sample$oxi_time), max(df_sample$oxi_time), length.out = 200)
  decay_curve <- data.frame(
    oxi_time = time_seq,
    CH4_pred = C0 * exp(-k * time_seq),
    sample_id = id
  )
  
  # Store
  curve_list[[id]] <- decay_curve
  data_list[[id]] <- df_sample
  meta_list[[id]] <- data.frame(sample_id = id, k = k, C0 = C0)
}

# Combine all
df_curve_all <- bind_rows(curve_list)
df_data_all <- bind_rows(data_list)
df_meta_all <- bind_rows(meta_list)

df_curve_all <- df_curve_all %>%
  left_join(df_data_all %>% select(sample_id, location) %>% distinct(), by = "sample_id")



# Initialize plot list
locations <- unique(df_data_all$location)

# Combine metadata into a dataframe
meta_list <- list()

for (id in unique(df_ch4$sample_id)) {
  df_sample <- df_ch4 %>% filter(sample_id == id)
  if (nrow(df_sample) < 3) next
  
  lm_fit <- lm(log(Amount) ~ oxi_time, data = df_sample)
  k <- -coef(lm_fit)[2]
  r2 <- summary(lm_fit)$r.squared
  C0 <- df_sample$Amount[1]
  
  # Store metadata
  meta_list[[id]] <- data.frame(sample_id = id, k = k, r2 = r2, C0 = C0)
}


# Create new labels: "10 (k = 0.0123)"
df_meta <- bind_rows(meta_list)

legend_labels <- setNames(
  paste0(" k = ", round(df_meta$k, 5),
         ", R² = ", round(df_meta$r2, 5)),
  df_meta$sample_id
)
my_pub_theme <- theme_classic() +
  theme(
    text = element_text(family = "Arial"),  # ← this line sets Arial globally
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(size = 5, margin = margin(t = 0.1, b=1)),
    axis.title.y = element_text(size = 5, margin = margin(r = 0.05)),
    axis.text = element_text(size = 5),
    
    # Force axis lines to be drawn
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    
    # Plot title & subtitle
    plot.subtitle = element_text(
      size = 6,
      #hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    plot.title = element_text(
      size = 5,
      face="bold",
      #  hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    
    # Legend text & title
    #  legend.title = element_text(size = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    
    # Legend spacing & key size
    legend.key.size = unit(0.1, "cm"),       # overall key box size
    # legend.key.height = unit(0.2, "cm"),     # vertical spacing
    #  legend.key.width = unit(0.2, "cm"),      # horizontal spacing
    legend.spacing.y = unit(0.07, "cm"),      # space between rows
    legend.spacing.x = unit(0.07, "cm"),      # space between columns
    legend.position = c(0.6, 0.9),
    legend.justification = c(0, 1),  
    # Legend background
    legend.background = element_rect(fill = "transparent")
  )

location_plots <- list()

for (loc in locations) {
  df_data_loc <- df_data_all %>% filter(location == loc)
  df_curve_loc <- df_curve_all %>% filter(location == loc)
  
  p <- ggplot() +
    geom_line(data = df_curve_loc, aes(x = oxi_time, y = CH4_pred, group = sample_id, color = sample_id),
              linewidth = 1, alpha=0.3) +
    geom_point(data = df_data_loc, aes(x = oxi_time, y = Amount, fill = sample_id),
               size = 1.5, color = "black", pch = 21, alpha = 1) +
    annotate("text",
             x = 260,  # right edge of plot
             y = 4,  # slightly above top data point
             label = "C(t) == C[0] %.% e^{-k * t}",
             parse = TRUE,
             hjust = 1.1, vjust = 2,
             size = 5 / .pt) +
    scale_fill_manual(values = custom_colors, labels = legend_labels)+ #, name = expression(C(t) == C[0] %.% e^{-k * t})) +
    scale_color_manual(values = custom_colors, labels = legend_labels )+ #, name = expression(C(t) == C[0] %.% e^{-k * t}))+
    labs(
      x = "Time [min]",
      y = expression(CH[4]~"[ppm]"),
      title = paste(loc)
    ) +
    my_pub_theme +
    coord_cartesian(
      xlim = c(0, 330),
      ylim = c(0, 4)
    )
  
  location_plots[[loc]] <- p
}

#location_plots[["Fussingø"]]
#wrap_plots(location_plots)
library(patchwork)
plot_comb<-location_plots[["Fussingø"]] + location_plots[["Silkeborg"]]

# ggsave("output/Fussingø_Silkeborg.svg",
#        plot_comb,
#        units = c("mm"),
#        height = 70,
#        width = 183,
#        dpi=300)
# 
# ggsave("output/Fussingø_Silkeborg_25_10_30.svg",
#        plot_comb,
#        units = c("mm"),
#        height = 60,
#        width = 160,
#        dpi=300)
# 


for (loc in locations) {
  df_data_loc <- df_data_all %>% filter(location == loc)
  df_curve_loc <- df_curve_all %>% filter(location == loc)
  
  p <- ggplot() +
    geom_line(data = df_curve_loc, aes(x = oxi_time, y = CH4_pred, group = sample_id, color = sample_id),
              linewidth = 1, alpha=0.3) +
    geom_point(data = df_data_loc, aes(x = oxi_time, y = Amount, fill = sample_id),
               size = 1.5, color = "black", pch = 21, alpha = 1) +
    annotate("text",
             x = 490,  # right edge of plot
             y = 3.1,  # slightly above top data point
             label = "C(t) == C[0] %.% e^{-k * t}",
             parse = TRUE,
             hjust = 1.1, vjust = 2,
             size = 5 / .pt) +
    scale_fill_manual(values = custom_colors, labels = legend_labels)+ #, name = expression(C(t) == C[0] %.% e^{-k * t})) +
    scale_color_manual(values = custom_colors, labels = legend_labels )+ #, name = expression(C(t) == C[0] %.% e^{-k * t}))+
    labs(
      x = "Time [min]",
      y = expression(CH[4]~"[ppm]"),
      title = paste(loc)
    ) +
    my_pub_theme +
    coord_cartesian(
      xlim = c(0, 600),
      ylim = c(0, 4)
    )+
    theme(legend.position = c(0.6, 0.7))
  
  location_plots[[loc]] <- p
}


plot_comb<-location_plots[["Dokkedal"]] 

# ggsave("output/Dokkedal.png",
#        plot_comb,
#        units = c("mm"),
#        height = 70,
#        width = 80,
#        dpi=300)
# 
# ggsave("output/Dokkedal.svg",
#        plot_comb,
#        units = c("mm"),
#        height = 70,
#        width = 80,
#        dpi=300)





###############################################################################
########################### Calculating specific rates ########################
###############################################################################

# Path to your Excel file
excel_file <- "data/CG_pressure_method.xlsx"
sheet_names <- excel_sheets(excel_file)

pressure_temp <- readxl::read_excel(excel_file, sheet = "Ark1",
                                    col_types = c("text", "numeric", "text", "date", "text", rep("guess", 20)))%>%
  slice(1:13) %>%      # keep rows 1 through 13
  select(1:5)          # keep the first 5 columns


# set the known volumes (in mL) for start and end headspace
V_start_mL <- 250   


# constants
R <- 8.314462618         # J / (mol K)
bar_to_Pa <- 1e5
mL_to_m3 <- 1e-6
m3_to_mL <- 1e6


pressure_temp_moles <- pressure_temp %>%
  mutate(
    Pressure_start = as.numeric(Pressure_start),
    Pressure_end = as.numeric(`Pressure at the end, assuming 2mL removal. Let sit at room temp.`), # adjust name if needed
    Temp_K = as.numeric(Temp) + 273.15,
    V_start = V_start_mL * mL_to_m3,
    V_end   = V_start,
    P_start_Pa = Pressure_start * bar_to_Pa,
    P_end_Pa   = Pressure_end * bar_to_Pa,
    n_start = ifelse(!is.na(P_start_Pa) & !is.na(Temp_K),
                     P_start_Pa * V_start / (R * Temp_K),
                     NA_real_),
    n_end = ifelse(!is.na(P_end_Pa) & !is.na(Temp_K),
                   P_end_Pa * V_end / (R * Temp_K),
                   NA_real_)
  ) %>%
  select(Sample, Pressure_start, Pressure_end, Date, Temp, n_start, n_end)%>%
  filter(!Sample=="neg")

df_ch4_pressure<-df_ch4%>%
  left_join(pressure_temp_moles, by =c("sample_id"="Sample"))%>%
  mutate(Temp=as.numeric(Temp))






# ---------------- Constants ----------------
V_headspace_mL <- 250    # headspace volume (mL)
V_sample_mL    <- 2      # sample withdrawn each time (mL)
#Temp_C_colname <- "Temp" # column name for temperature in °C
Pressure_colname <- "Pressure_start" # column name with measured headspace pressure (bar)
Pressure_end_colname <- "Pressure_end" # column name with end pressure (bar) to anchor backward recon
# ------------------------------------------------

# constants / conversions
R <- 8.314462618 # J/mol*K
bar_to_Pa <- 1e5
mL_to_m3 <- 1e-6
m3_to_mL <- 1e6

# derived volumes in m3
V_head_m3 <- V_headspace_mL * mL_to_m3
V_sample_m3 <- V_sample_mL * mL_to_m3

# functions to convert between n and P (P in bar, V in m3, T in K)
PVT_to_n <- function(P_bar, V_m3, Temp_K) { (P_bar * bar_to_Pa) * V_m3 / (R * Temp_K) }

nVT_to_Pbar <- function(n_mol, V_m3, Temp_K) { (n_mol * R * Temp_K / V_m3) / bar_to_Pa }

# Prepare grouped data for reconstruction: compute Temp_K, index per injection,
# and compute measured n at start/end when pressures are available
dat_grouped <- df_ch4_pressure %>%
  mutate(Temp_K = Temp + 273.15) %>%
  group_by(sample_id) %>%
  arrange(oxi_time, `Injection Name`, .by_group = TRUE) %>%
  mutate(idx = row_number(),                            
    n_start_measured = if_else(!is.na(Pressure_start) & !is.na(Temp_K),
                         PVT_to_n(Pressure_start, V_head_m3, Temp_K),
                         NA_real_),
    # measured moles at end pressure (for backward anchor) per row
    n_end_measured = if_else(!is.na(Pressure_end) & !is.na(Temp_K),
                             PVT_to_n(Pressure_end, V_head_m3, Temp_K),
                             NA_real_)
  ) %>%
  ungroup()

# Process per sample_id: reconstruct forward and backward gas volumes and pressures
results <- dat_grouped %>%
  group_by(sample_id) %>%
  group_modify(~ {
    d <- .x
    n <- nrow(d)
    if (n == 0) return(d)  # skip empty groups
    
    # Initialize result vectors
    n_forward      <- rep(NA_real_, n)
    P_forward_bar  <- rep(NA_real_, n)
    n_back         <- rep(NA_real_, n)
    P_back_bar     <- rep(NA_real_, n)
    
    # ================================
    # --- Forward reconstruction ---
    # ================================
    # Start from measured initial gas amount (n_start_measured)
    if (!is.na(d$n_start_measured[1])) {
      n_forward[1]     <- d$n_start_measured[1]
      P_forward_bar[1] <- d$Pressure_start[1]
      
      # Iterate forward over rows
      if (n >= 2) {
        for (i in 2:n) {
          if (!is.na(n_forward[i - 1])) {
            # Calculate gas removed when taking a 2 mL sample at previous pressure
            n_removed <- PVT_to_n(P_forward_bar[i - 1], V_sample_m3, d$Temp_K[i - 1])
            
            # Remaining moles of gas after removal
            n_after <- n_forward[i - 1] - n_removed
            
            # Store updated values
            n_forward[i]     <- n_after
            P_forward_bar[i] <- nVT_to_Pbar(n_after, V_head_m3, d$Temp_K[i])
          }
        }
      }
    }
    
    # ================================
    # --- Backward reconstruction ---
    # ================================
    # Start from measured final gas amount (n_end_measured)
    if (!is.na(d$n_end_measured[n])) {
      n_back[n]     <- d$n_end_measured[n]
      P_back_bar[n] <- d$Pressure_end[n]
      
      # Iterate backward over rows
      if (n >= 2) {
        for (i in (n - 1):1) {
          if (!is.na(n_back[i + 1])) {
            # Determine pressure to use for withdrawn sample
            P_withdraw_bar <- d$Pressure_start[i]
            
            # Fallbacks: try next start or final end pressure if missing
            if (is.na(P_withdraw_bar)) {
              P_withdraw_bar <- d$Pressure_start[i + 1]
              if (is.na(P_withdraw_bar)) P_withdraw_bar <- d$Pressure_end[n]
            }
            
            if (!is.na(P_withdraw_bar)) {
              # Convert withdrawn gas to moles
              n_sample <- PVT_to_n(P_withdraw_bar, V_sample_m3, d$Temp_K[i])
              
              # Gas amount before removal
              n_before <- n_back[i + 1] + n_sample
              
              # Store updated values
              n_back[i]     <- n_before
              P_back_bar[i] <- nVT_to_Pbar(n_before, V_head_m3, d$Temp_K[i])
            } else {
              # If pressure still unknown, carry NA forward
              n_back[i]     <- NA_real_
              P_back_bar[i] <- NA_real_}}}}}
    
    # =======================================
    # Return augmented data frame for group
    # =======================================
    d %>%
      mutate(
        n_forward              = n_forward,
        P_forward_bar          = P_forward_bar,
        n_removed_forward      = lag(n_forward) - n_forward,
        V_removed_mL_forward   = V_sample_mL,
        n_back                 = n_back,
        P_back_bar             = P_back_bar,
        n_added_backward       = lead(n_back) - n_back,
        V_sample_mL_report     = V_sample_mL
      )
  }) %>%
  ungroup()


# Keep only columns you want to inspect and ensure order is preserved
final <- results %>%
  arrange(sample_id, idx) %>%
  select(
    sample_id, idx, `Injection Name`, oxi_time,
    Pressure_start, Pressure_end, Temp_K,
    n_start_measured, n_forward, P_forward_bar, n_removed_forward, V_removed_mL_forward,
    n_end_measured, n_back, P_back_bar, n_added_backward, V_sample_mL_report
  )

print(final, n = Inf)



# 1) Row-wise differences added to final
final_comp <- final %>%
  mutate(
    delta_n = n_forward - n_back,                 # positive => forward predicts more moles
    delta_P_bar = P_forward_bar - P_back_bar,     # positive => forward predicts higher pressure
    rel_delta_n = if_else(!is.na(n_back) & n_back != 0, delta_n / n_back, NA_real_),
    rel_delta_P = if_else(!is.na(P_back_bar) & P_back_bar != 0, delta_P_bar / P_back_bar, NA_real_)
  )

# 2) Per-sample summary statistics
summary_by_sample <- final_comp %>%
  group_by(sample_id) %>%
  summarise(
    n_rows = n(),
    mean_delta_n = mean(delta_n, na.rm = TRUE),
    sd_delta_n = sd(delta_n, na.rm = TRUE),
    max_abs_delta_n = max(abs(delta_n), na.rm = TRUE),
    rmse_delta_n = sqrt(mean((delta_n)^2, na.rm = TRUE)),
    mean_delta_P_bar = mean(delta_P_bar, na.rm = TRUE),
    max_abs_delta_P_bar = max(abs(delta_P_bar), na.rm = TRUE),
    final_delta_n = last(delta_n),                # difference at the last index in group
    final_delta_P_bar = last(delta_P_bar)
  ) %>%
  arrange(desc(max_abs_delta_n))

# 3) Quick global diagnostics
global_diag <- final_comp %>%
  summarise(
    total_rows = n(),
    rows_with_both = sum(!is.na(n_forward) & !is.na(n_back)),
    overall_rmse_n = sqrt(mean((delta_n)^2, na.rm = TRUE)),
    overall_mean_abs_delta_n = mean(abs(delta_n), na.rm = TRUE),
    overall_mean_rel_delta_n = mean(abs(rel_delta_n), na.rm = TRUE)
  )

# 4) Show worst mismatches (top 10 by absolute mole difference)
worst_rows <- final_comp %>%
  arrange(desc(abs(delta_n))) %>%
  slice_head(n = 10) %>%
  select(sample_id, idx, `Injection Name`, oxi_time,
         Pressure_bar, Pressure_end_bar, Temp_K,
         n_measured, n_forward, n_back, delta_n, delta_P_bar, V_sample_mL_report)

# Print results
print(global_diag)
print(summary_by_sample, n = Inf)
print(worst_rows, n = Inf)

### So, the max difference between begin and end is 0.018 bar. This is quite ok I think.
# Combine forward/backward reconstructions: use the mean when both estimates exist,
# otherwise use the single available estimate (used to compute CH4 moles below).
final_combined <- final %>%
  mutate(
    # combined mole count: mean when both present, otherwise the one that exists
    n_combined = case_when(
      !is.na(n_forward) & !is.na(n_back) ~ (n_forward + n_back) / 2,
      !is.na(n_forward) & is.na(n_back)  ~ n_forward,
      is.na(n_forward) & !is.na(n_back)  ~ n_back,
      TRUE                                ~ NA_real_
    ),
    
    # optional: combined pressure computed from n_combined (bar)
    P_combined_bar = if_else(
      !is.na(n_combined) & !is.na(Temp_K),
      (n_combined * R * Temp_K / V_head_m3) / bar_to_Pa,
      NA_real_
    ),
    
    # useful flag to know how the value was formed
    n_source = case_when(
      !is.na(n_forward) & !is.na(n_back) ~ "mean(forward,back)",
      !is.na(n_forward) & is.na(n_back)  ~ "forward_only",
      is.na(n_forward) & !is.na(n_back)  ~ "backward_only",
      TRUE                                ~ "none"
    )
  )

# inspect
print(final_combined %>% select(sample_id, idx, `Injection Name`, oxi_time,
                                n_forward, n_back, n_combined, P_combined_bar, n_source),
      n = Inf)



# Compute CH4 molecule counts per injection using combined headspace n, convert to µmol,
# and compute natural log for linear fitting
df_ch4_moles <- df_ch4 %>%
  left_join(final_combined, by = c("sample_id", "oxi_time", "Injection Name")) %>%
  mutate(
    n_ch4 = (Amount / 1e6) * n_combined, 
    n_ch4_umol = n_ch4 * 1e6,  
    ln_n_ch4_umol = log(n_ch4_umol)
  )


# Constants
R <- 8.314462618         # J/mol/K
V_head_m3 <- 250 * 1e-6  # 250 mL → m³
Temp_atm_K <- 15 + 273.15    # 15°C
Pressure_bar_atm <- 1.013     # bar absolute
atm_ch4 <- 1.93           # ppm  https://gml.noaa.gov/ccgg/trends_ch4/

PVT_to_n <- function(P_bar, V_m3, Temp_K) { (P_bar * bar_to_Pa) * V_m3 / (R * Temp_K) }

C0_mol_ch4 <- (atm_ch4/1e6)*PVT_to_n(Pressure_bar_atm, V_head_m3, Temp_atm_K)
C0_µmol_ch4 <-C0_mol_ch4*1e6 


# Per-sample log-linear fit of ln(n_CH4[µmol]) versus time to extract first-order rate k
df_slopes <- df_ch4_moles %>%
  group_by(sample_id) %>%
  summarise(
    # Fit ln(n_ch4) ~ oxi_time
    slope = {
      fit <- lm(ln_n_ch4_umol ~ oxi_time)
      coef(fit)[2]   # extract slope
    },
    intercept = {
      fit <- lm(ln_n_ch4_umol ~ oxi_time)
      coef(fit)[1]   # extract intercept
    },
    r2 = summary(lm(ln_n_ch4_umol ~ oxi_time))$r.squared,
    m_soil = first(m_soil),
    .groups = "drop"
  ) %>%
  mutate(k = -slope) %>%  # since ln(n) = a - k*t
  mutate(k_hour=k*60)%>%
  mutate(r=k_hour*C0_µmol_ch4)%>%
  mutate(r_soil=r/m_soil)



############# Now to plotting ################
# Build plotmath legend labels per sample_id using k, r_soil and R^2 (used in ggplot legends)
legend_labels_rates <- setNames(
  lapply(1:nrow(df_slopes), function(i) {
    k_val <- round(df_slopes$k_hour[i], 3)
    r_val <- df_slopes$r_soil[i]
    r2_val <- round(df_slopes$r2[i], 4)
    
    # Convert r_val to scientific notation with 2 decimals
    r_val_sci <- formatC(r_val, format = "e", digits = 2)
    
    # Extract mantissa and exponent for plotmath
    parts <- strsplit(r_val_sci, "e")[[1]]
    mantissa <- as.numeric(parts[1])
    exponent <- as.integer(parts[2])
    
    # Build plotmath expression
    bquote("k =" ~ .(k_val) * " " * h^-1 * ", " *
             r[soil] == .(mantissa) %.% 10^.(exponent) * " " * mu * mol ~ h^-1 ~ g^-1 * ", " *
             R^2 == .(r2_val))
  }),
  df_slopes$sample_id
)


my_pub_theme <- theme_classic() +
  theme(
    text = element_text(family = "Arial"),  # ← this line sets Arial globally
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(size = 5, margin = margin(t = 0.1, b=1)),
    axis.title.y = element_text(size = 5, margin = margin(r = 0.05)),
    axis.text = element_text(size = 5),
    
    # Force axis lines to be drawn
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    
    # Plot title & subtitle
    plot.subtitle = element_text(
      size = 6,
      #hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    plot.title = element_text(
      size = 5,
      face="bold",
      #  hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    
    # Legend text & title
    #  legend.title = element_text(size = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    
    # Legend spacing & key size
    legend.key.size = unit(0.1, "cm"),       # overall key box size
    # legend.key.height = unit(0.2, "cm"),     # vertical spacing
    #  legend.key.width = unit(0.2, "cm"),      # horizontal spacing
    legend.spacing.y = unit(0.07, "cm"),      # space between rows
    legend.spacing.x = unit(0.07, "cm"),      # space between columns
    legend.position = c(0.35, 1),
    legend.justification = c(0, 1),  
    # Legend background
    legend.background = element_rect(fill = "transparent")
  )



df_ch4_neg <- df_all %>%
  filter(Component == "CH4 [FID]") %>%
#  filter(!is.na(location))%>%
  filter(oxi_time < 600) %>%
  arrange(sample_id, oxi_time)%>%
  filter(sample_id=="neg")

location_plots <- list()

for (loc in locations) {
  df_data_loc <- df_data_all %>% filter(location == loc)
  df_curve_loc <- df_curve_all %>% filter(location == loc)
  
  p <- ggplot() +
    # negative control line (dashed black, no grouping)
    geom_line(data = df_ch4_neg, 
      aes(x = oxi_time, y = Amount),        # change these column names if needed
      linewidth = 0.8,
      inherit.aes = FALSE, 
      color = "grey80") +
    geom_point(
      data = df_ch4_neg, 
      aes(x = oxi_time, y = Amount),        # change these column names if needed
      fill = "grey80",size = 1.5, color = "black", pch = 21) +
    geom_line(data = df_curve_loc, aes(x = oxi_time, y = CH4_pred, group = sample_id, color = sample_id),
              linewidth = 1, alpha=0.3) +
    geom_point(data = df_data_loc, aes(x = oxi_time, y = Amount, fill = sample_id),
               size = 1.5, color = "black", pch = 21, alpha = 1) +
    annotate("text",
             x = 180,  # right edge of plot
             y = 4.42,  # slightly above top data point
             label = "C(t) == C[0] %.% e^{-k * t}",
             parse = TRUE,
             hjust = 1.1, vjust = 2,
             size = 5 / .pt) +
    scale_fill_manual(values = custom_colors, labels = legend_labels_rates)+ #, name = expression(C(t) == C[0] %.% e^{-k * t})) +
    scale_color_manual(values = custom_colors, labels = legend_labels_rates )+ #, name = expression(C(t) == C[0] %.% e^{-k * t}))+
    labs(
      x = "Time [min]",
      y = expression(CH[4]~"[ppm]"),
      title = paste(loc)
    ) +
    my_pub_theme +
    coord_cartesian(
      xlim = c(0, 400),
      ylim = c(0, 4)
    )
  
  location_plots[[loc]] <- p
}



for (loc in locations) {
  df_data_loc <- df_data_all %>% filter(location == loc)
  df_curve_loc <- df_curve_all %>% filter(location == loc)
  
  p <- ggplot() +
    geom_line(data = df_curve_loc, aes(x = oxi_time, y = CH4_pred, group = sample_id, color = sample_id),
              linewidth = 1, alpha=0.3) +
    geom_point(data = df_data_loc, aes(x = oxi_time, y = Amount, fill = sample_id),
               size = 1.5, color = "black", pch = 21, alpha = 1) +
    annotate("text",
             x = 330,  # right edge of plot
             y = 1.1,  # slightly above top data point
             label = "C(t) == C[0] %.% e^{-k * t}",
             parse = TRUE,
             hjust = 1.1, vjust = 2,
             size = 5 / .pt) +
    scale_fill_manual(values = custom_colors, labels = legend_labels_rates)+ #, name = expression(C(t) == C[0] %.% e^{-k * t})) +
    scale_color_manual(values = custom_colors, labels = legend_labels_rates )+ #, name = expression(C(t) == C[0] %.% e^{-k * t}))+
    labs(
      x = "Time [min]",
      y = expression(CH[4]~"[ppm]"),
      title = paste(loc)
    ) +
    my_pub_theme +
    coord_cartesian(
      xlim = c(0, 600),
      ylim = c(0, 4)
    )+
    theme(legend.position = c(0.35, 0.25))
  
  location_plots[[loc]] <- p
}




plot_comb<-location_plots[["Dokkedal"]] 

ggsave("oxi_rates/Dokkedal_rates_15C.png",
       plot_comb,
       units = c("mm"),
       height = 70,
       width = 80,
       dpi=300)
# 
ggsave("oxi_rates/Dokkedal_rates.svg",
       plot_comb,
       units = c("mm"),
       height = 70,
       width = 80,
       dpi=300)





################## Combining them #######################
# Subset curves/data for the two locations to be combined in a single panel
df_curve_capsa<-df_curve_all%>%filter(location %in% c("Fussingø", "Silkeborg"))
df_all_capsa<-df_data_all%>%filter(location %in% c("Fussingø", "Silkeborg"))

df_ch4_neg <- df_all %>%
  filter(Component == "CH4 [FID]") %>%
  filter(oxi_time < 600) %>%
  arrange(sample_id, oxi_time)%>%
  filter(sample_id=="neg")





# 1) Get final CH4_pred per sample_id (largest oxi_time)
df_final <- df_curve_capsa %>%
  group_by(sample_id, location) %>%
  slice_max(order_by = oxi_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(location, sample_id, oxi_time, CH4_pred) %>%
  rename(final_oxi_time = oxi_time, final_CH4_pred = CH4_pred)

# 2) Build mapping ordered by descending final_CH4_pred within each location
mapping <- df_final %>%
  group_by(location) %>%
  arrange(location, desc(final_CH4_pred)) %>%   # highest CH4_pred first
  mutate(
    seq = row_number(),
    loc_letter = toupper(substr(location, 1, 1)),
    sample_id_v2 = paste0(loc_letter, seq)
  ) %>%
  ungroup() %>%
  select(location, sample_id, sample_id_v2, final_CH4_pred)

# 3) Join mapping back to your full dataset (so every row gets sample_id_v2)
df_curve_capsa2 <- df_curve_capsa %>%
  left_join(mapping %>% select(location, sample_id, sample_id_v2), by = c("location", "sample_id"))

# 4) If you want only the endpoint rows with the new labels:
df_curve_labels <- df_curve_capsa2 %>%
  group_by(sample_id) %>%
  slice_max(order_by = oxi_time, n = 1, with_ties = FALSE) %>%
  ungroup()


df_curve_labels_add<-df_curve_labels%>%select(sample_id, sample_id_v2)



# negative control last point
df_ctrl_label <- df_ch4_neg %>%
  slice_max(order_by = oxi_time, n = 1, with_ties = FALSE) %>%
  mutate(sample_id = "Neg control")   # label text

# Join labels for remapped sample names (e.g., F1, F2, S1...) used for plotting
df_slopes2<-df_slopes%>%
  left_join(df_curve_labels)%>%
  filter(!is.na(sample_id_v2))

#### Making the legends:
legend_labels_rates <- setNames(
  lapply(seq_len(nrow(df_slopes2)), function(i) {
    # values
    sample_label <- df_slopes2$sample_id_v2[i]
    k_val <- round(df_slopes2$k_hour[i], 3)
    r_val <- df_slopes2$r_soil[i]
    r2_val <- round(df_slopes2$r2[i], 4)
    
    # scientific notation with 2 decimals
    r_val_sci <- formatC(r_val, format = "e", digits = 2)
    parts <- strsplit(r_val_sci, "e")[[1]]
    mantissa <- as.numeric(parts[1])
    exponent <- as.integer(parts[2])
    
    # build plotmath expression; sample_label inserted as plain text
    bquote(.(sample_label) * ":" ~
             k == .(k_val) * " " * h^-1 * "," ~
             r[soil] == .(mantissa) %.% 10^.(exponent) * " " * mu * mol ~ h^-1 ~ g^-1 * "," ~
             R^2 == .(r2_val))
  }),
  df_slopes2$sample_id_v2
)




df_curve_capsa<-df_curve_capsa%>%left_join(df_curve_labels_add)
df_all_capsa<-df_all_capsa%>%left_join(df_curve_labels_add)


custom_colors <- c(
  "F1" =   "brown", # f7081d
  "F2" =  "brown", # Purple
  "F3" =  "brown", # Coral
  "F4" =   "brown",   # gold
  
  "S1" ="#2219e6", # Teal
  "S2" = "#2219e6", # Periwinkle
  "S3" ="#2219e6", # Blue
  "S4" ="#2219e6"
)


legend_title_expr <- expression(C(t) == C[0] %.% e^{-k * t})


p <- ggplot() +
    # negative control line (dashed black, no grouping)
    geom_line(data = df_ch4_neg, 
              aes(x = oxi_time, y = Amount),        # change these column names if needed
              linewidth = 0.8,
              inherit.aes = FALSE, 
              color = "grey80") +
    geom_point(
      data = df_ch4_neg, 
      aes(x = oxi_time, y = Amount),        # change these column names if needed
      fill = "grey80",size = 1.3, color = "black", pch = 21) +
    geom_line(data = df_curve_capsa, aes(x = oxi_time, y = CH4_pred, group = sample_id_v2, color = sample_id_v2),
              linewidth = 1, alpha=0.3) +
    geom_point(data = df_all_capsa, aes(x = oxi_time, y = Amount, fill = sample_id_v2),
               size = 1.3, color = "black", pch = 21, alpha = 1) +
  geom_text(data = df_curve_labels,
                  aes(x = oxi_time, y = CH4_pred, label = sample_id_v2),
                  #nudge_x = 8,            # push label to the right
                  #direction = "y",
                  hjust = -0.5,
               #   vjust=0.5,
                  #segment.size = 0.2,
                  size = 5 / .pt,
                  show.legend = FALSE) +
  # label for negative control
  geom_text(data = df_ctrl_label,
                  aes(x = oxi_time, y = Amount, label = sample_id),
                  color = "grey40",
                  hjust = -0.1,
                  size = 5 / .pt,
                  show.legend = FALSE) +
    scale_fill_manual(name = legend_title_expr, values = custom_colors, labels = legend_labels_rates)+ #, name = expression(C(t) == C[0] %.% e^{-k * t})) +
    scale_color_manual(name = legend_title_expr, values = custom_colors, labels = legend_labels_rates )+ #, name = expression(C(t) == C[0] %.% e^{-k * t}))+
    labs(
      x = "Time [min]",
      y = expression(CH[4]~"[ppm]")
    ) +
  guides(
    color = guide_legend(ncol = 1, byrow = F),
    fill  = guide_legend(ncol = 1, byrow = F))+
    my_pub_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 5))+
    coord_cartesian(
      xlim = c(0, 420),
      ylim = c(0, 4),  clip = "off"
    )

p

# 
# ggsave("output/combined_oxi.png",
#        p,
#        units = c("mm"),
#        height = 80,
#        width = 180,
#        dpi=300)
# # 
# ggsave("output/combined_oxi.svg",
#        p,
#        units = c("mm"),
#        height = 80,
#        width = 180,
#        dpi=300)
# # 




####### Adding table next to plot ##############
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
tbl <- df_slopes2 %>%
  dplyr::select(sample_id_v2, k_hour, r_soil, r2) %>%
  dplyr::mutate(
    Sample = sample_id_v2,
    k = sprintf("%.2f", round(k_hour, 3)),
    r_soil = formatC(r_soil, format = "e", digits = 3),
    R2 = sprintf("%.3f", round(r2, 3))
  ) %>%
  dplyr::select(Sample, k, r_soil, R2) %>%
  arrange(Sample)

# create tableGrob with plain headers (we will replace them)
tbl_theme <- ttheme_minimal(
  core = list(fg_params = list(fontsize = 5)),      # 5 pt body
  colhead = list(fg_params = list(fontface = "bold", fontsize = 5)) # initial header style
)

tbl_grob <- tableGrob(tbl, rows = NULL, theme = tbl_theme)

# define the plotmath expressions for the headers
# Define header expressions with units using plotmath
# - k [h^-1]
# - r[soil] [10^-4 * mu * mol ~ h^-1 ~ g^-1]
header_exprs <- list(
  expression("Sample"),
  expression(k ~ "[" ~ h^-1 ~ "]"),
  expression(r[soil] ~ "[" ~ mu * mol ~ h^-1 ~ g^-1 ~ "]"),
  expression(R^2)
)


# find the indices of the column header grobs in the tableGrob
# the layout names for header text are typically "colhead-fg"
colhead_idx <- which(tbl_grob$layout$name == "colhead-fg")

# sanity check: if not found, print available layout names (helps debugging)
if (length(colhead_idx) != length(header_exprs)) {
  message("Unexpected header grob layout. Available layout names:")
  print(unique(tbl_grob$layout$name))
}

# replace each header grob with a textGrob using the expression
for (i in seq_along(header_exprs)) {
  if (i <= length(colhead_idx)) {
    tg <- textGrob(label = header_exprs[[i]],
                   gp = gpar(fontface = "bold", fontsize = 5),
                   x = unit(0.5, "npc"), just = "centre")
    tbl_grob$grobs[[colhead_idx[i]]] <- tg
  }
}

# optional: tweak column widths
tbl_grob$widths <- unit(c(0.4, 0.4, 0.7, 0.4), "in")
#tbl_grob$heights <- tbl_grob$heights * 0.8

# remove legend from plot and combine
p_no_legend <- p + theme(legend.position = "none")

final <- plot_grid(
  p_no_legend,
  ggdraw(tbl_grob),
  ncol = 2,
  rel_widths = c(3, 1)
)

print(final)

ggsave("output/combined_oxi_15C.png",
       final,
       units = c("mm"),
       height = 80,
       width = 183,
       dpi=300)
# 
ggsave("output/combined_oxi_15C_26_01_20.svg",
       final,
       units = c("mm"),
       height = 80,
       width = 183,
       dpi=300)
# 



############################# less wide no legend ##########################


ggsave("output/combined_oxi_25_12_08.svg",
       p_no_legend,
       units = c("mm"),
       height = 60,
       width = 80,
       dpi=300)
