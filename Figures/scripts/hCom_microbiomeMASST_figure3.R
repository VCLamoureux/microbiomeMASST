setwd("/Users/vincentlamoureux/Library/CloudStorage/")

library(readxl)
library(tidyverse)
library(pheatmap)
library(mixOmics)
library(vegan)
library(caret)
library(ggpubr)
library(stringr)
library(ComplexUpset)
library(duckplyr)

# import tables
feature_table <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/mzmine/fbmn_V2_quant.csv") |>
  dplyr::rename_with(~gsub(" Peak area$", "", .x)) |>
  dplyr::select(-tidyselect::matches("^Blank_(7|8|9)$")) |>
  dplyr::rename(`three_alpha_12_K.mzML` = `three_alpha_12-K.mzML`, Feature = `row ID`) |>
  dplyr::mutate(Feature = as.character(Feature))

# import quality control table
drift <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/mzmine/RT_hCom_V2.csv")
colnames(drift)[1] <- "Feature"
metadata <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/metadata_hCom_1_15052024.csv")
annotation <- readr::read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/annotations_hCom_V2.tsv") |>
  dplyr::rename(Feature = `#Scan#`) |>
  dplyr::mutate(Feature = as.character(Feature))

# import the randomized sequence
sequence <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/randomized_Sequence_hCom1_13052024.csv")
colnames(sequence) <- sequence[1, ]
sequence <- sequence[-1, ]
sequence_1 <- sequence |> 
  dplyr::mutate(order = dplyr::row_number(),
                Time  = dplyr::case_when(str_detect(`Sample ID`, "Plate") ~ 48,
                                         str_detect(`Sample ID`, "regrow") ~ 72,
                                         TRUE ~ NA_real_)) |>
  dplyr::rename(filename = `File Name`) |>
  dplyr::select(order, Time, tidyselect::everything())

# get info feature
info_feature_name <- feature_table |>
  dplyr::select(Feature, mz = 2, RT = 3) |>
  dplyr::left_join(annotation, by = "Feature")

# transpose data
data_transpose <- feature_table |>
  tibble::column_to_rownames("Feature") |>
  dplyr::select(ends_with(".mzML")) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("filename")

# assess retention time drift 
rt_data_ordered <- drift |>
  dplyr::filter(Feature == 23462) |>
  dplyr::select(contains("Feature RT")) |>
  tidyr::pivot_longer(everything(), names_to = "filename", values_to = "RT") |>
  dplyr::mutate(RT = as.numeric(RT),
                filename = sub("\\.mzML Feature RT$", "", filename)) |>
  dplyr::filter(filename != "0") |>
  dplyr::left_join(sequence_1, by = "filename") |>
  dplyr::filter(!str_detect(filename, "ACN|Blank|three|TCA|TDCA|TLCA|TCDCA")) |>
  dplyr::arrange(order)

rt_stats <- rt_data_ordered |>
  dplyr::summarise(Mean_RT = mean(RT, na.rm = TRUE),
                   SD_RT   = stats::sd(RT, na.rm = TRUE),
                   Var_RT  = stats::var(RT, na.rm = TRUE))

ggplot(rt_data_ordered, aes(x = stats::reorder(filename, order), y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "Retention time (min)", title = "Retention Time Drift (Feature 23462)")

# evaluate peak intensity drift
peakarea_data_ordered <- drift |>
  dplyr::filter(Feature == 23462) |>
  dplyr::select(contains("Peak")) |>
  tidyr::pivot_longer(everything(), names_to = "filename", values_to = "Peak") |>
  dplyr::mutate(filename = sub("\\.mzML Peak area$", "", filename)) |>
  dplyr::left_join(sequence_1, by = "filename") |>
  dplyr::filter(!str_detect(filename, "ACN|Blank|three|TCA|TDCA|TLCA|TCDCA|sixmix")) |>
  dplyr::arrange(order)

peakarea_stats <- peakarea_data_ordered |>
  dplyr::summarise(Mean_peakarea = mean(Peak, na.rm = TRUE),
                   SD_peakarea   = stats::sd(Peak, na.rm = TRUE),
                   Var_peakarea  = stats::var(Peak, na.rm = TRUE))

ggplot(peakarea_data_ordered, aes(x = stats::reorder(filename, order), y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "Peak area", title = "Peak Area Drift (Feature 23462)")

# Check QCmix samples
data_sixmix <- data_transpose |>
  dplyr::filter(str_detect(filename, "sixmix")) |>
  dplyr::mutate(filename = sub("\\.mzML$", "", filename), filename = sub(" Peak area$", "", filename))

sixmix_feature_info <- tibble::tibble(
  Feature     = colnames(data_sixmix)[-1],
  Mean_sixmix = data_sixmix |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_sixmix   = data_sixmix |> tibble::column_to_rownames("filename") |> apply(2, stats::sd)) |>
  dplyr::mutate(CV_sixmix = SD_sixmix / Mean_sixmix) |>
  dplyr::filter(Mean_sixmix > 0) |>
  dplyr::arrange(dplyr::desc(Mean_sixmix))

sixmix_feature_info_name <- sixmix_feature_info |>
  dplyr::left_join(info_feature_name, by = "Feature") |>
  dplyr::select(Feature, Mean_sixmix, SD_sixmix, CV_sixmix, mz, RT, Compound_Name)

# Check QCpool samples
data_pool <- data_transpose |>
  dplyr::filter(str_detect(filename, "pool")) |>
  dplyr::mutate(filename = sub("\\.mzML$", "", filename))

pool_feature_info <- tibble::tibble(
  Feature   = colnames(data_pool)[-1],
  Mean_pool = data_pool |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_pool   = data_pool |> tibble::column_to_rownames("filename") |> apply(2, stats::sd)) |>
  dplyr::mutate(CV_pool = SD_pool / Mean_pool) |>
  dplyr::filter(Mean_pool > 0) |>
  dplyr::arrange(dplyr::desc(Mean_pool))

pool_feature_info_name <- pool_feature_info |>
  dplyr::left_join(info_feature_name, by = "Feature")

# Check Blank samples
data_blank <- data_transpose |>
  dplyr::filter(str_detect(filename, "Blank"))

blank_feature_info <- tibble::tibble(
  Feature    = colnames(data_blank)[-1],
  Mean_blank = data_blank |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_blank   = data_blank |> tibble::column_to_rownames("filename") |> apply(2, stats::sd)) |>
  dplyr::mutate(CV_blank = SD_blank / Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |>
  dplyr::arrange(dplyr::desc(Mean_blank))

blank_feature_info_name <- blank_feature_info |>
  dplyr::left_join(info_feature_name, by = "Feature")

# check the samples
data_sample <- data_transpose |>
  dplyr::filter(!str_detect(filename, "ACN|TCDCA|TCA|TDCA|TLCA|three_alpha|sixmix|pool|Blank")) |>
  dplyr::mutate(filename = paste0("sample_", filename))

sample_feature_info <- tibble::tibble(
  Feature     = colnames(data_sample)[-1],
  Mean_sample = data_sample |> tibble::column_to_rownames("filename") |> colMeans(),
  SD_sample   = data_sample |> tibble::column_to_rownames("filename") |> apply(2, stats::sd)) |>
  dplyr::mutate(CV_sample = SD_sample / Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |>
  dplyr::arrange(dplyr::desc(Mean_sample))

sample_feature_info_name <- sample_feature_info |>
  dplyr::left_join(info_feature_name, by = "Feature")

# remove specific features matching the 6 stds using in the QCmix
specified_features <- c("14800", "23462", "15635", "17817", "28028", "30197")

feature_to_remove <- blank_feature_info |>
  dplyr::left_join(pool_feature_info, by = "Feature") |>
  dplyr::filter(Mean_blank > 0) |>
  dplyr::mutate(pool_Blank = Mean_pool / Mean_blank) |>
  dplyr::filter(pool_Blank < 5 | is.na(pool_Blank)) |>
  dplyr::bind_rows(blank_feature_info |> 
                     dplyr::filter(Feature %in% specified_features)) |>
  dplyr::distinct(Feature, .keep_all = TRUE)

data_clean <- data_transpose |>
  dplyr::select(-any_of(feature_to_remove$Feature))

# remove features that are not 5 times higher in pool than QCmix
feature_to_remove_qcmix <- sixmix_feature_info |>
  dplyr::left_join(pool_feature_info, by = "Feature") |>
  dplyr::filter(Mean_sixmix > 0) |>
  dplyr::mutate(pool_Mix = Mean_pool / Mean_sixmix) |>
  dplyr::filter(pool_Mix < 5 | is.na(pool_Mix)) |>
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

data_clean2 <- data_clean |>
  dplyr::select(-any_of(feature_to_remove_qcmix$Feature)) |>
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank|ACN"))

# transpose back the df
data_clean_transpose <- data_clean2 |>
  tibble::column_to_rownames("filename") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Feature")

# Keep only bile acid features matching to known and unknown bile acids libraries 
features_to_include <- annotation |>
  dplyr::filter(LibraryName %in% c("BILELIB19.mgf", "GNPS-BILE-ACID-MODIFICATIONS.mgf")) |>
  dplyr::pull(Feature) |>
  unique()

# merge with annotation
merged <- dplyr::left_join(data_clean_transpose, annotation, by = "Feature")

# pivot to long format
df_long <- merged |>
  dplyr::select(Feature, tidyselect::ends_with(".mzML")) |>
  dplyr::filter(Feature %in% features_to_include) |>
  tidyr::pivot_longer(cols = -Feature, names_to = "sample", values_to = "intensity") |>
  dplyr::mutate(sample     = sub("\\.mzML$", "", sample),
                condition  = sub("_(\\d+)$", "", sample),
                replicate  = sub("^.*_(\\d+)$", "\\1", sample),
                replicate  = factor(replicate),
                Feature    = as.character(Feature))

# combine control samples
df_long_std <- df_long |>
  dplyr::mutate(condition_std = dplyr::case_when(condition %in% c("BHI", "media_mPYG_BHI") ~ "BHI", TRUE ~ condition))

# filter for C. spiroforme and A. hallii vs BHI
df_cs <- df_long_std |>
  dplyr::filter(condition_std %in% c("Clostridium_spiroforme_DSM_1552", "BHI"))
df_eh <- df_long_std |>
  dplyr::filter(condition_std %in% c("Eubacterium_hallii_DSM_3353", "BHI"))

# keep features that are at least 3 times higher in C. spiroforme than in BHI
kept_features <- df_cs |>
  dplyr::group_by(Feature, condition_std) |>
  dplyr::summarise(m = mean(intensity, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = condition_std, values_from = m, values_fill = 0, names_repair = "minimal") |>
  dplyr::filter(`Clostridium_spiroforme_DSM_1552` >= 3 * `BHI`) |>
  dplyr::pull(Feature)

# keep only relevant columns from annotation
annotation_small <- annotation |>
  dplyr::select(Feature, LibraryName, Compound_Name, SpecMZ, SpectrumID)

# replace df_cs for df_eh for the other microbes and Eubacterium_hallii_DSM_3353
df_cs_1 <- df_eh |>
  dplyr::filter(Feature %in% kept_features, condition_std == "Eubacterium_hallii_DSM_3353") |>
  dplyr::left_join(annotation_small, by = "Feature")

# calculate the mean of the replicates
df_cs_1_mean <- df_cs_1 |>
  dplyr::group_by(Feature, condition_std) |>
  dplyr::summarise(mean_intensity = mean(intensity, na.rm = TRUE), .groups = "drop_last") |>
  dplyr::ungroup()

# join with annotation
df_increase_cs <- df_cs_1_mean |>
  dplyr::group_by(Feature) |>
  dplyr::left_join(annotation_small, by = "Feature")

# remove water losses adducts for visualization, create delta mass and hydroxylation columns
df_cs_1_mean_final <- df_cs_1_mean |>
  dplyr::semi_join(df_increase_cs, by = "Feature") |>
  dplyr::left_join(annotation_small, by = "Feature") |>
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = dplyr::if_else(str_detect(Compound_Name, "delta mass"),
                                            str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"), Compound_Name)) |>
  dplyr::mutate(hydroxylation = dplyr::case_when(
    str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
    str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
    str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
    str_detect(Compound_Name, regex("dihydroxy", ignore_case = TRUE)) ~ "di",
    str_detect(Compound_Name, regex("aMCA|bMCA", ignore_case = TRUE)) ~ "tri",
    str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
    str_detect(Compound_Name, regex("12-hydroxy|3-hydroxy", ignore_case = TRUE)) ~ "mono",
    str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
    str_detect(Compound_Name, regex("deoxycholic", ignore_case = TRUE)) ~ "di",
    str_detect(Compound_Name, regex("litho", ignore_case = TRUE)) ~ "mono",
    str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
      !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
    str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~
      tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[, 2]),
    TRUE ~ NA_character_)) |>
  dplyr::filter(mean_intensity > 0)

cat_levels <- c("mono","di","tri","tetra","penta","unknown")

# split bile acid libraries
df_cs_bile <- df_cs_1_mean_final |>
  dplyr::filter(LibraryName == "BILELIB19.mgf", mean_intensity > 0)
df_cs_gnps <- df_cs_1_mean_final |>
  dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf", mean_intensity > 0)

# metadata map: delta_mass -> hydroxylation
comp_meta <- dplyr::bind_rows(
  df_cs_bile |> dplyr::select(delta_mass, hydroxylation),
  df_cs_gnps |> dplyr::select(delta_mass, hydroxylation)) |>
  dplyr::distinct() |>
  dplyr::mutate(hydroxylation = tidyr::replace_na(hydroxylation, "unknown"))

hmap <- stats::setNames(comp_meta$hydroxylation, comp_meta$delta_mass)

# order function: by class (cat_levels), then alphabetical
order_by_class <- function(x) {
  cls <- hmap[x]
  cls[is.na(cls)] <- "unknown"
  primary  <- match(cls, cat_levels)
  primary[is.na(primary)] <- length(cat_levels)
  secondary <- rank(tolower(x), ties.method = "first")
  x[order(primary, secondary)]
}

bilelib_compounds <- unique(df_cs_bile$delta_mass)
gnps_compounds    <- unique(df_cs_gnps$delta_mass)

bile_order <- order_by_class(bilelib_compounds)
gnps_order <- order_by_class(gnps_compounds)

center_node <- "Clostridium spiroforme"
group_bile  <- "BILELIB19"
group_gnps  <- "GNPS"

edges_center <- tibble::tibble(from = c(center_node, center_node),
                               to   = c(group_bile, group_gnps))

edges_bile   <- tibble::tibble(from = rep(group_bile, length(bile_order)),
                               to   = bile_order)

edges_gnps   <- tibble::tibble(from = rep(group_gnps, length(gnps_order)),
                               to   = gnps_order)

edges        <- dplyr::bind_rows(edges_center, edges_bile, edges_gnps)

g <- igraph::graph_from_data_frame(edges, directed = TRUE)

hyd_vec <- unname(hmap[igraph::V(g)$name])
hyd_vec[is.na(hyd_vec)] <- "unknown"
igraph::V(g)$hydroxylation <- hyd_vec

pal_h <- c(mono="#96D5F7", di="#1EC6BD", tri="#FFE070", tetra="#FF9F45", penta="#FF6F91", unknown="grey80")

p <- ggraph::ggraph(g, layout = "dendrogram", circular = TRUE) +
  ggraph::geom_edge_diagonal(color = "gray40", alpha = 0.6) +
  ggraph::geom_node_point(ggplot2::aes(filter = (name == center_node)),
                          size = 10, color = "darkred", stroke = NA) +
  ggraph::geom_node_point(ggplot2::aes(filter = (name %in% c(group_bile, group_gnps))),
                          size = 8, color = "darkorange", stroke = NA) +
  ggraph::geom_node_point(ggplot2::aes(filter = (!name %in% c(center_node, group_bile, group_gnps)),
                                       color  = hydroxylation),
                          size = 5.5, stroke = NA) +
  ggplot2::scale_color_manual(values = pal_h, na.value = "grey80") +
  ggraph::geom_node_text(
    ggplot2::aes(label  = name,
                 color  = hydroxylation,
                 filter = leaf & !(name %in% c(center_node, group_bile, group_gnps)),
                 angle  = -ggraph::node_angle(x, y)),
    hjust = "outward", size = 2.6, show.legend = FALSE
  ) +
  ggplot2::guides(color = "none") +
  ggplot2::theme_void() +
  ggplot2::ggtitle("Radial Network: Clostridium spiroforme → BILELIB19 & GNPS (grouped by hydroxylation)")

p

#ggsave("Radial_Network_Ah_hydroxylation.pdf", plot = p, width = 11, height = 11, dpi = 900) 
#getwd()

p <- ggraph::ggraph(g, layout = "dendrogram", circular = TRUE) +
  ggraph::geom_edge_diagonal(color = "gray40", alpha = 0.6) +
  ggraph::geom_node_point(aes(filter = (name == center_node)), size = 10, color = "darkred", stroke = NA) +
  ggraph::geom_node_point(aes(filter = (name %in% c(group_bile, group_gnps))), size = 8, color = "darkorange", stroke = NA) +
  ggraph::geom_node_point(aes(filter = (!name %in% c(center_node, group_bile, group_gnps)), color = hydroxylation),
                          size = 5.5, stroke = NA) +
  ggplot2::scale_color_manual(values = pal_h, na.value = "grey80") +
  ggplot2::guides(color = "none") +
  ggplot2::theme_void() +
  ggplot2::ggtitle("Radial Network: Clostridium spiroforme → BILELIB19 & GNPS (grouped by hydroxylation)")
p
#ggsave("Radial_Network_Ah_hydroxylation_1.pdf", plot = p, width = 11, height = 11, dpi = 900) 
#getwd()

