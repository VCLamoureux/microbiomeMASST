# set working directory
setwd("/Users/vincentlamoureux/Downloads/")

# import libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(janitor)
library(purrr)
library(tibble)
library(patchwork)


file <- "ace_PILOT_DORRESTEON (Modified)_20251030_171226.csv"
enc <- readr::guess_encoding(file, n_max = 2000)$encoding[1]
raw_lines <- readr::read_lines(file, locale = readr::locale(encoding = enc))
header_row <- which(str_detect(raw_lines, "\\bA1\\b\\s*,\\s*Time \\[ms\\]"))[1]

df <- readr::read_csv(file, skip = header_row - 1, locale = readr::locale(encoding = enc), name_repair = "minimal") |> 
  clean_names()

nm <- names(df)
time_idx <- which(str_detect(nm, "^time_ms"))
signal_idx <- time_idx - 1

pairs <- tibble(signal = nm[signal_idx], time = nm[time_idx])

#create the tibble
tidy_dat <- map2_dfr(pairs$signal, pairs$time, ~{
  tibble(
    well = .x,
    time_min = as.numeric(df[[.y]]) / 60000,
    value = as.numeric(df[[.x]]))}) |> 
  dplyr::filter(!is.na(time_min), !is.na(value))

dose_map <- tribble(
  ~well, ~conc_nM, ~series,
  
  "b3", 0.1,  "ENP",
  "b4", 0.2,  "ENP",
  "b5", 0.5,  "ENP",
  "b6", 1.0,  "ENP",
  "b7", 1.5,  "ENP",
  "b8", 2.0,  "ENP",
  "b9", 2.5,  "ENP",
  "b10", 3.0, "ENP",
  #"b11", 3.5, "ENP",
  "c3", 5.0,  "ENP",
  "c4", 10.0, "ENP",
  "c5", 20.0, "ENP",
  
  "d3", 0.1,  "NP_ENP",
  "d4", 0.2,  "NP_ENP",
  "d5", 0.5,  "NP_ENP",
  "d6", 1.0,  "NP_ENP",
  "d7", 1.5,  "NP_ENP",
  "d8", 2.0,  "NP_ENP",
  "d9", 2.5,  "NP_ENP",
  "d10", 3.0, "NP_ENP",
  #"d11", 3.5, "NP_ENP",
  "e3", 5.0,  "NP_ENP",
  "e4", 10.0, "NP_ENP",
  "e5", 20.0, "NP_ENP",
  
  "f3", NA, "EC",
  "f4", NA, "EC",
  "f5", NA, "EC",
  
  "f6", 0.1, "POS_IC",
  "f7", 0.1, "POS_IC",
  "f8", 0.1, "POS_IC",
  
  "g3", NA, "BC",
  "g5", NA, "BC",
  "g6", NA, "BC",
  "g7", NA, "BC")

# Join to tidy dataset
tidy_dat_1 <- tidy_dat |> left_join(dose_map, by = "well")

# controls EC for enzyme control and BC for background control
EC_wells <- c("f3","f4","f5")
BC_wells <- c("g3","g5","g6","g7")

# Subtract background
bc_trace <- tidy_dat_1 |>
  dplyr::filter(well %in% BC_wells) |>
  dplyr::group_by(time_min) |>
  dplyr::summarise(bc_mean = mean(value, na.rm = TRUE), .groups = "drop")

# Continuous background function; rule=2 clamps outside observed range
bc_fun <- with(bc_trace, approxfun(time_min, bc_mean, rule = 2))

tidy_dat_bc <- tidy_dat_1 |>
  dplyr::mutate(value_bc = value - bc_fun(time_min))

# Calculate slope
slope_10_30 <- function(df){
  d <- df |>
    dplyr::filter(time_min >= 10, time_min <= 30) |>
    dplyr::arrange(time_min)
  if (nrow(d) < 3) return(tibble(slope = NA_real_, r2 = NA_real_))
  fit <- lm(value_bc ~ time_min, data = d)
  tibble(
    slope = coef(fit)[["time_min"]],
    r2    = round(summary(fit)$r.squared * 100, 1)
  )
}

well_slopes <- tidy_dat_bc |>
  dplyr::group_by(well, series, conc_nM) |>
  dplyr::group_modify(~slope_10_30(.x)) |>
  dplyr::ungroup()

# EC normalization & % inhibition 
EC_mean_slope <- well_slopes |>
  dplyr::filter(well %in% EC_wells) |>
  dplyr::summarise(m = mean(slope, na.rm = TRUE)) |> dplyr::pull(m)

# clamp the value between 0 and 100
clamp01 <- function(x) pmin(pmax(x, 0), 100)

results_well <- well_slopes |>
  dplyr::mutate(
    pct_inhibition_raw = 100 * (EC_mean_slope - slope) / EC_mean_slope,
    pct_inhibition     = clamp01(pct_inhibition_raw),
    pct_activity       = 100 - pct_inhibition)


# plot data
plot_df <- results_well |>
  dplyr::filter(!is.na(series), !is.na(conc_nM), !series %in% c("POS_IC")) |>
  dplyr::mutate(conc_lab = factor(conc_nM, levels = sort(unique(conc_nM))))
#write_csv(plot_df, "enalaprilat_enzyme_assay_results_final.csv")

plot_df_enalaprilat <- ggplot(plot_df, aes(x = conc_lab, y = pct_inhibition, color = series)) +
  geom_line(aes(group = series), linewidth = 1) +
  geom_point(size = 3.5, shape = 19, alpha = 0.9, stroke = NA) +
  scale_color_manual(values = c(ENP = "#00aeef", NP_ENP = "#ec008c")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Concentration (nM)", y = "% Relative Inhibition", color = "Series") +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 16), axis.title = element_text(face = "bold"), legend.position = "right")

plot_df_enalaprilat
#ggsave("enalaprilat_enzyme_assay.pdf", plot = plot_df_enalaprilat, width = 6, height = 5, dpi = 900, bg = "transparent")
#getwd()

# choose which signal to show on traces:
use_value <- "value"

trace_dat <- tidy_dat_1 |>
  dplyr::filter(series %in% c("ENP","NP_ENP"), conc_nM %in% c(0.5, 5, 10, 20)) |>
  dplyr::mutate(conc_lab = factor(conc_nM, levels = c(0.5, 5, 10, 20)))

# Endpoints
endpts <- trace_dat |>
  dplyr::group_by(series, conc_lab, well) |>
  dplyr::slice_max(time_min, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::group_by(series, conc_lab) |>
  dplyr::summarise(
    time_min = max(time_min, na.rm = TRUE),
    y        = median(.data[[use_value]], na.rm = TRUE),
    .groups  = "drop") |> dplyr::mutate(lbl = paste0(as.character(conc_lab), " nM"))

# create color palettes
enp_palette    <- colorRampPalette(c("#a6e7ff", "#00aeef"))(4)
npenp_palette  <- colorRampPalette(c("#f8b7dc", "#ec008c"))(4)

shade_cols <- c(
  "ENP_0.5" = enp_palette[1],
  "ENP_5"   = enp_palette[2],
  "ENP_10"  = enp_palette[3],
  "ENP_20"  = enp_palette[4],
  "NP_ENP_0.5" = npenp_palette[1],
  "NP_ENP_5"   = npenp_palette[2],
  "NP_ENP_10"  = npenp_palette[3],
  "NP_ENP_20"  = npenp_palette[4])

trace_dat <- trace_dat |> dplyr::mutate(series_conc = paste(series, conc_lab, sep = "_"))
#write_csv(trace_dat, "enalaprilat_enzyme_assay_traces.csv")

p_traces <- ggplot(
  trace_dat, aes(x = time_min,
      y = .data[[use_value]],
      color = series_conc,
      group = interaction(series_conc, well))) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = shade_cols, name = "Group × Concentration") +
  geom_text(data = endpts,
    aes(x = time_min + 0.5, y = y, label = lbl, color = paste(series, conc_lab, sep = "_")), inherit.aes = FALSE, size = 4, fontface = "bold", hjust = 0) +
  coord_cartesian(xlim = c(0, 60), clip = "off") +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  labs(x = "Time (minutes)", y = "Absorbance (OD)") +
  theme_classic(base_size = 18) +
  theme(plot.margin = margin(10, 60, 10, 10), plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold"), legend.position = "none")

p_traces
#ggsave("enalaprilat_gradient_traces.pdf",plot = p_traces, width = 5, height = 5, dpi = 900, bg = "transparent")

# add a single value
v30 <- trace_dat |>
  dplyr::group_by(series, conc_lab) |>
  dplyr::summarise( y30 = approx(x = time_min, y = value, xout = 30, rule = 2, ties = mean)$y, .groups = "drop") |>
  dplyr::mutate(series_conc = paste(series, conc_lab, sep = "_"), lbl30 = sprintf("%.3f", y30))

ggplot(v30, aes(x = conc_lab, y = y30, color = series)) +
  geom_point(position = position_dodge(width = 0.35), size = 3.5, shape = 19) +
  scale_color_manual(values = c(ENP = "#00aeef", NP_ENP = "#ec008c")) +
  scale_y_continuous(limits = c(0.7, 1.0), breaks = seq(0.7, 1.0, by = 0.1)) +
  labs(
    x = "Concentration (nM)",
    y = "Absorbance at 30 min (OD)",
    color = "Series") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"))
#ggsave("enalaprilat_30min_values.pdf", width = 4, height = 5, dpi = 900, bg = "transparent") 












