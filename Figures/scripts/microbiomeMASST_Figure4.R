# Set working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# load packages 
library(tidyverse)
library(reshape2)

# remove mzml and peak area extensions
strip_mzml <- function(x) gsub("\\.mzML$", "", x)
rm_peak_area <- function(x) gsub(" Peak area", "", x)

# files to excluse as they have been rerun
file_to_exclude <- c(
  "group3_B2T0.mzML","group3_B3T72.mzML","group3_B6T0.mzML","group3_B7T72.mzML",
  "group3_D10T0.mzML","group3_D11T72.mzML","group3_D2T0.mzML",
  "group3_D5T0.mzML","group3_D6T0.mzML","group3_D9T72.mzML",
  "pool_qc_groups_26.mzML","pool_qc_groups_27.mzML","sixmix_24.mzML")

exclude_stds <- c("GCDCA_20uM.mzML","GCA_20uM.mzML","TDCA_20uM.mzML", "TCA_20uM.mzML","CA_20uM.mzML")

# import metadata
metadata <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/Drug_screening/Drug_metabolites/metadata_drug_metabolites.csv")
metadata <- metadata |>
  dplyr::mutate( ATTRIBUTE_time   = as.numeric(stringr::str_extract(filename, "(?<=T)\\d+")), ATTRIBUTE_growth = paste0(ATTRIBUTE_time, "_", ATTRIBUTE_groups)) |>
  dplyr::filter(!filename %in% EXCLUDE_FILES) |>
  dplyr::mutate(
    ATTRIBUTE_drug_mixture = dplyr::case_when(
      stringr::str_detect(filename, "B2|C2|D2|E2|F2|G2")   ~ "drug_mix_1",
      stringr::str_detect(filename, "B3|C3|D3|E3|F3|G3")   ~ "drug_mix_2",
      stringr::str_detect(filename, "B4|C4|D4|E4|F4|G4")   ~ "drug_mix_3",
      stringr::str_detect(filename, "B5|C5|D5|E5|F5|G5")   ~ "drug_mix_4",
      stringr::str_detect(filename, "B6|C6|D6|E6|F6|G6")   ~ "drug_mix_5",
      stringr::str_detect(filename, "B7|C7|D7|E7|F7|G7")   ~ "drug_mix_6",
      stringr::str_detect(filename, "B8|C8|D8|E8|F8|G8")   ~ "drug_mix_7",
      stringr::str_detect(filename, "B9|C9|D9|E9|F9|G9")   ~ "drug_mix_8",
      stringr::str_detect(filename, "B10|C10|D10|E10|F10|G10") ~ "drug_mix_9",
      stringr::str_detect(filename, "B11|C11|D11|E11|F11|G11") ~ "drug_mix_10",
      TRUE ~ ATTRIBUTE_drug_mixture))

# import randomized sequence
randomized_sequence <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/Drug_screening/Drug_metabolites/normal_flow_randomized_10_10_2024_order.csv")
colnames(randomized_sequence) <- as.character(unlist(randomized_sequence[1, ]))
randomized_sequence <- randomized_sequence[-1, ]
colnames(randomized_sequence)[2] <- "filename"
randomized_sequence$Order <- as.numeric(randomized_sequence$Order)

# import quality control
drift <- readr::read_csv(
  "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/Drug_screening/Drug_metabolites/Normal_flow/MZmine_output/output/others_quant_mzmine")
colnames(drift)[1] <- "Feature"

feature_table <- readr::read_csv(
  "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/Drug_screening/Drug_metabolites/Normal_flow/MZmine_output/output/others_iimn_fbmn_quant_drug_metabolites_2.csv") |>
  dplyr::rename_all(rm_peak_area)

annotations <- readr::read_tsv(
  "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/Drug_screening/Drug_metabolites/Normal_flow/bdd2bea069514fa7871cfe517c388df3/nf_output/library/merged_results_with_gnps.tsv")
annotations$`#Scan#` <- as.character(annotations$`#Scan#`)
annotations <- annotations |> dplyr::rename(Feature = `#Scan#`)

info_feature <- feature_table[, 1:3]
colnames(info_feature) <- c("Feature", "mz", "RT")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_name <- dplyr::left_join(info_feature, annotations, by = "Feature")

# transpose table
data_transpose <- feature_table |>
  tibble::column_to_rownames("row ID") |>
  dplyr::select(dplyr::contains(".mzML")) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("filename") |>
  dplyr::filter(!filename %in% file_to_exclude, !filename %in% exclude_stds)

# check the RT drift of the internal standard
feature_row <- drift |> dplyr::filter(Feature == 43281)
rt_columns  <- feature_row |> dplyr::select(dplyr::contains("RT"))

rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT")
rt_data$filename <- gsub(".mzML Feature RT", "", rt_data$filename)

rt_data_ordered <- rt_data |>
  dplyr::left_join(randomized_sequence, by = "filename") |>
  dplyr::filter(!filename %in% c(
    "GCDCA_20uM","GCA_20uM","TDCA_20uM","TCA_20uM","CA_20uM",
    "Blank_1","pool_qc_groups_26_20241026094150")) |>
  dplyr::mutate(Order = as.numeric(Order)) |>
  dplyr::arrange(Order)

rt_stats <- rt_data |>
  dplyr::summarise(
    Mean_RT     = mean(RT, na.rm = TRUE),
    SD_RT       = sd(RT, na.rm = TRUE),
    Variance_RT = var(RT, na.rm = TRUE))

ggplot(rt_data_ordered, aes(x = reorder(filename, Order), y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 9183")

# check the intensity drift of the internal standard
feature_row_peakarea <- drift |> dplyr::filter(Feature == 43281)
peakarea_columns <- feature_row_peakarea |> dplyr::select(dplyr::contains("Peak"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data$filename <- gsub(".mzML Peak area", "", peakarea_data$filename)

peakarea_data_ordered <- peakarea_data |>
  dplyr::left_join(randomized_sequence, by = "filename") |>
  dplyr::arrange(Order) |>
  dplyr::filter(!stringr::str_detect(filename, "pool|sixmix"))

peakarea_stats <- peakarea_data |>
  dplyr::summarise(
    Mean_peakarea = mean(Peak, na.rm = TRUE),
    SD_peakarea   = sd(Peak, na.rm = TRUE),
    Variance_peakarea = var(Peak, na.rm = TRUE))

ggplot(peakarea_data_ordered, aes(x = reorder(filename, Order), y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "peakarea", title = "peakarea Drift for Feature 9183")

# check the QCmix samples
data_sixmix <- data_transpose |> dplyr::filter(stringr::str_detect(filename, "sixmix"))
data_sixmix$filename <- strip_mzml(data_sixmix$filename)
data_sixmix$filename <- gsub("Peak area", "", data_sixmix$filename)

sixmix_feature_info <- data.frame(
  Feature     = colnames(data_sixmix)[-1],
  Mean_sixmix = data_sixmix |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_sixmix   = data_sixmix |> tibble::column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sixmix = SD_sixmix / Mean_sixmix) |>
  dplyr::filter(Mean_sixmix > 0) |>
  dplyr::arrange(dplyr::desc(Mean_sixmix))

sixmix_feature_info_name <- dplyr::left_join(
  sixmix_feature_info, info_feature_name, by = "Feature") |>
  dplyr::select("Feature","Mean_sixmix","SD_sixmix","CV_sixmix","mz","RT","Compound_Name")

# Pool
data_pool <- data_transpose |> dplyr::filter(stringr::str_detect(filename, "pool"))
data_pool$filename <- strip_mzml(data_pool$filename)
data_pool$filename <- gsub("Peak area", "", data_pool$filename)

pool_feature_info <- data.frame(
  Feature   = colnames(data_pool)[-1],
  Mean_pool = data_pool |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_pool   = data_pool |> tibble::column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_pool = SD_pool / Mean_pool) |>
  dplyr::filter(Mean_pool > 0) |>
  dplyr::arrange(dplyr::desc(Mean_pool))

# Blank
data_blank <- data_transpose |> dplyr::filter(stringr::str_detect(filename, "Blank"))
blank_feature_info <- data.frame(
  Feature    = colnames(data_blank)[-1],
  Mean_blank = data_blank |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_blank   = data_blank |> tibble::column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank / Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |>
  dplyr::arrange(dplyr::desc(Mean_blank))

# Samples only
data_sample <- data_transpose |>
  dplyr::left_join(metadata |> dplyr::select(filename, ATTRIBUTE_type), by = "filename") |>
  dplyr::filter(ATTRIBUTE_type %in% c("sample","Control")) |>
  dplyr::select(-ATTRIBUTE_type)

sample_feature_info <- data.frame(
  Feature      = colnames(data_sample)[-1],
  Mean_sample  = data_sample |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_sample    = data_sample |> tibble::column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample / Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |>
  dplyr::arrange(dplyr::desc(Mean_sample))

pool_feature_info_name   <- dplyr::left_join(pool_feature_info,   info_feature_name, by = "Feature")
sample_feature_info_name <- dplyr::left_join(sample_feature_info, info_feature_name, by = "Feature")
blank_feature_info_name  <- dplyr::left_join(blank_feature_info,  info_feature_name, by = "Feature")

# remove features 
feature_to_remove <- blank_feature_info |>  
  left_join(pool_feature_info) |> 
  dplyr::filter(Mean_blank > 0) |> 
  dplyr::mutate(pool_Blank = Mean_pool/Mean_blank) |> 
  dplyr::filter(pool_Blank < 5 | is.na(pool_Blank)) |> 
  dplyr::bind_rows(blank_feature_info |> dplyr::filter(Feature %in% specified_features)) |> 
  dplyr::distinct(Feature, .keep_all = TRUE)

# Remove any features listed in 'feature_to_remove'
data_clean <- data_transpose |> 
  dplyr::select(-c(feature_to_remove$Feature))

# Features to be removed samples/QCmix < 5
feature_to_remove_qcmix <- sixmix_feature_info |> 
  left_join(pool_feature_info) |> 
  dplyr::filter(Mean_sixmix > 0) |> 
  dplyr::mutate(pool_Mix = Mean_pool/Mean_sixmix) |> 
  dplyr::filter(pool_Mix < 5 | is.na(pool_Mix)) |> 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal, pool, blanks
data_clean2 <- data_clean |> 
  dplyr::select(-c(feature_to_remove_qcmix$Feature)) |> 
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank"))

data_clean_transpose <- data_clean2 |>
  tibble::column_to_rownames("filename") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Feature")

# pivot longer
df_melted <- reshape2::melt(
  data_clean_transpose,
  id.vars = "Feature",
  variable.name = "filename",
  value.name = "Peak_Area")
df_melted$filename <- strip_mzml(df_melted$filename)
df_melted_filtered <- df_melted |>
  dplyr::filter(!stringr::str_detect(filename, "sixmix|Blank|pool")) |>
  dplyr::mutate(
    Group           = stringr::str_extract(filename, "^[^_]+"),
    ATTRIBUTE_time  = as.numeric(stringr::str_extract(filename, "(?<=T)\\d+")),
    ATTRIBUTE_growth = paste0(Group, "_", ATTRIBUTE_time))

metadata_filtered <- metadata
metadata_filtered$filename <- strip_mzml(metadata_filtered$filename)

df_melted_plot <- df_melted_filtered |>
  dplyr::filter(Feature == "37650") |>
  dplyr::left_join(metadata_filtered |> dplyr::select(filename, ATTRIBUTE_drug_mixture), by = "filename") |>
  dplyr::filter(ATTRIBUTE_drug_mixture %in% c("drug_mix_6"))

# Calculate the t-test
tt_res <- df_melted_plot %>%
  dplyr::mutate(
    Peak_Area      = as.numeric(Peak_Area),
    ATTRIBUTE_time = as.numeric(ATTRIBUTE_time)) %>%
  dplyr::filter(ATTRIBUTE_time %in% c(0, 72), !is.na(Peak_Area)) %>%
  dplyr::group_by(Feature, Group) %>%
  dplyr::group_modify(~{
    dat <- .x
    n0  <- sum(dat$ATTRIBUTE_time == 0,  na.rm = TRUE)
    n72 <- sum(dat$ATTRIBUTE_time == 72, na.rm = TRUE)
    m0  <- mean(dat$Peak_Area[dat$ATTRIBUTE_time == 0],  na.rm = TRUE)
    m72 <- mean(dat$Peak_Area[dat$ATTRIBUTE_time == 72], na.rm = TRUE)
    
    if (n0 >= 2 && n72 >= 2) {
      tt <- t.test(Peak_Area ~ factor(ATTRIBUTE_time), data = dat, var.equal = TRUE)
      tibble::tibble(
        n_0 = n0, n_72 = n72,
        mean_0 = m0, mean_72 = m72,
        diff_72_minus_0 = m72 - m0,
        t_stat = unname(tt$statistic),
        df     = unname(tt$parameter),
        p_value = tt$p.value
      )
    } else {
      tibble::tibble(
        n_0 = n0, n_72 = n72,
        mean_0 = m0, mean_72 = m72,
        diff_72_minus_0 = m72 - m0,
        t_stat = NA_real_, df = NA_real_, p_value = NA_real_
      )
    }
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(Feature, Group)

tt_res

# colors to plot
custom_colors        <- c("#4DB3E6","#FFB74D","#A5D6A7","#FF8A65","#9575CD")
custom_colors_darker <- c("#1976D2","#F57C00","#388E3C","#D84315","#5E35B1")

df_melted_plot_1 <- ggplot(df_melted_plot, aes(x = ATTRIBUTE_growth, y = Peak_Area, fill = Group)) +
  geom_boxplot(aes(fill = Group), color = "black", size = 0.5, alpha = 0.5, width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color = Group), shape = 19, position = position_jitter(0.15), size = 6, alpha = 1, stroke = NA) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors_darker) +
  labs(x = "Condition", y = "Raw peak area") +
  theme_minimal(base_size = 22) +
  theme(
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x  = element_text(angle = 90, vjust = 1, hjust = 1, size = 22),
    axis.text.y  = element_text(size = 22),
    legend.position = "none",
    strip.text   = element_text(size = 8),
    axis.line    = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks   = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"))

df_melted_plot_1
#ggsave("enalapril_proline.pdf", df_melted_plot_1, width = 5, height = 7, dpi = 300)










