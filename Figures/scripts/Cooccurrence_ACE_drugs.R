# set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth")

library(tidyverse)
library(networkD3)
library(duckplyr)
library(htmlwidgets)
library(ComplexUpset)

#  point towards all the generated files
files <- list.files(
  path       = "/Users/vincentlamoureux/Desktop/enalapril_microbiomeMASST_FASST_2",
  pattern    = "\\.csv$",
  full.names = TRUE)

# load all the files
all_matches <- files |>
  purrr::map_dfr(~{
    compound_name <- basename(.x) |>
      str_remove(".csv$") |>
      str_remove("_")
    read_csv(.x) |>
      dplyr::mutate(compound = compound_name)
  })

# formatting and filtering
all_matches_datasets <- all_matches |> 
  dplyr::mutate(
    dataset        = sub("^[^:]+:([^:]+):.*$", "\\1", USI),
    filename_level = str_remove(USI, ":scan:.*$")) |>
  dplyr::filter(`Matching Peaks` >= 4, Cosine >= 0.7)

# check all molecules have been loaded at least 1 time
all_matches_datasets_unique <- all_matches_datasets |>
  dplyr::distinct(compound_name, .keep_all = TRUE)

# remove all non-human datasets
df <- all_matches_datasets |>
  dplyr::filter(!dataset %in% c(
    "MSV000094528", "MSV000093469", "MSV000083593", "MSV000083772", "MSV000084908",
    "MSV000085843", "MSV000087612", "MSV000087685", "MSV000092754", "MSV000094642",
    "NORMAN-61a72f15-591b-5cf7-9a4e-7e26549f1564", "MSV000092487", "MSV000089190",
    "MSV000096359", "MSV000097575", "MSV000097976", "NORMAN-b3198128-7d7f-5d05-86a1-60da4b920140", "ST001715", "MSV000096860", "MSV000083889"))

all_compounds <- sort(unique(df$compound))

# match compound names to be used in the sankey
enalaprilat         <- "Enalaprilat"
enalaprilat_proline <- "enalaprilat_proline"
enalapril_proline   <- "enalapril_proline"

# layer 4 final nodes
final_nodes <- c(enalaprilat_proline, "no_proline", "no_proline_3")

# prodrugs (every compound except Enalaprilat, enalaprilat_proline, enalapril_proline)
prodrugs <- setdiff(all_compounds, c(enalaprilat, enalaprilat_proline, enalapril_proline))

# check for enalaprilat and desprolyl_enalaprilat
file_flags <- df |>
  group_by(filename_level) |>
  summarise(
    has_enalaprilat = any(compound == enalaprilat),
    has_proline     = any(compound == enalaprilat_proline),
    .groups         = "drop")

# filenames that contain Enalaprilat
enalap_files <- file_flags |>
  filter(has_enalaprilat) |>
  pull(filename_level)

# counts for Enalaprilat + proline vs Enalaprilat no proline
n_enalap_and_proline <- file_flags |>
  filter(has_enalaprilat, has_proline) |>
  nrow()

n_enalap_no_proline <- file_flags |>
  filter(has_enalaprilat, !has_proline) |>
  nrow()

# build links (edges) for the Sankey plot
# dataset -> prodrug
parent_metabolite <- df |>
  dplyr::filter(compound %in% prodrugs) |>
  distinct(dataset, filename_level, compound) |>
  count(dataset, compound, name = "value") |>
  dplyr::rename(source = dataset, target = compound)

# prodrug -> Enalaprilat using ALL files with Enalaprilat
prodrug_to_enalaprilat <- df |>
  dplyr::filter(filename_level %in% enalap_files, compound %in% prodrugs) |>
  distinct(filename_level, compound) |>
  count(compound, name = "value") |>
  dplyr::rename(source = compound) |>
  dplyr::mutate(target = enalaprilat)

# which prodrugs actually co-occur with Enalaprilat
prodrugs_with_enalaprilat <- unique(prodrug_to_enalaprilat$source)

# prodrugs that do not co-occur with Enalaprilat... create dummy node no_proline_2
prodrugs_no_enalaprilat <- setdiff(prodrugs, prodrugs_with_enalaprilat)

dummy_prodrugs <- tibble(
  source = prodrugs_no_enalaprilat,
  target = "no_proline_2",
  value  = 1)

# no_proline_2 -> no_proline_3
link_no_proline2_to_3 <- tibble(
  source = "no_proline_2",
  target = "no_proline_3",
  value  = sum(dummy_prodrugs$value))

# Enalaprilat -> enalaprilat_proline
link_enalaprilat_to_proline <- tibble(
  source = enalaprilat,
  target = enalaprilat_proline,
  value  = n_enalap_and_proline)

# Enalaprilat -> no_proline
link_enalaprilat_to_no_proline <- tibble(
  source = enalaprilat,
  target = "no_proline",
  value  = n_enalap_no_proline)

# enalapril (2nd layer) -> enalapril_proline (3rd layer)
# count filenames where both enalapril and enalapril_proline appear
enalapril_proline_files <- df |>
  dplyr::filter(compound == enalapril_proline) |>
  pull(filename_level) |>
  unique()

enalapril_to_enalapril_proline <- df |>
  dplyr::filter(filename_level %in% enalapril_proline_files, compound == "enalapril") |>
  distinct(filename_level) |>
  summarise(value = n()) |>
  dplyr::mutate(source = "enalapril", target = enalapril_proline)

# combine all links
links <- bind_rows(
  parent_metabolite,
  prodrug_to_enalaprilat,      
  dummy_prodrugs,             
  link_no_proline2_to_3,         
  link_enalaprilat_to_proline,   
  link_enalaprilat_to_no_proline,
  enalapril_to_enalapril_proline )

# Build nodes (layers) - only ones that actually appear in links
dataset_names <- sort(unique(parent_metabolite$source))

# prodrugs that actually appear in edges
prodrugs_used <- sort(
  unique(c(parent_metabolite$target,
           prodrug_to_enalaprilat$source,
           dummy_prodrugs$source,
           enalapril_to_enalapril_proline$source)))

nodes_ordered <- c(dataset_names, prodrugs_used, c(enalaprilat, enalapril_proline, "no_proline_2"), final_nodes)

nodes <- data.frame(name = nodes_ordered, stringsAsFactors = FALSE)

# map node names to indices for networkD3
links <- links |>
  dplyr::mutate(
    IDsource = match(source, nodes$name) - 1,
    IDtarget = match(target, nodes$name) - 1)

# Generate Sankey plot
sn <- sankeyNetwork(
  Links       = links,
  Nodes       = nodes,
  Source      = "IDsource",
  Target      = "IDtarget",
  Value       = "value",
  NodeID      = "name",
  sinksRight  = TRUE,
  iterations  = 0,
  nodeWidth   = 35,
  nodePadding = 10,
  fontSize    = 14,
  width       = 500,
  height      = 1000)

sn

html_file <- "Enalaprilat_sankey.html"
#saveWidget(sn, file = html_file, selfcontained = TRUE)
#browseURL(html_file)


# Calculate odds ratio and p-value and 95 CI 
# Parent drugs 
ace_drugs <- c("enalapril", "Enalaprilat","ramipril", "quinapril", "moexipril", "trandolapril")

# drug metabolites
despro_compounds <- c("enalapril_proline", "enalaprilat_proline")

# Filter the df for these compounds
ace_df <- df %>% dplyr::filter(compound %in% c(ace_drugs, despro_compounds))

# generate a true and false based on the presense of enalapril(at), desprolyl compounds, and the other ACEs
ace_flags <- ace_df %>%
  group_by(filename_level) %>%
  summarise(
    has_enalapril_group = any(compound %in% c("enalapril", "Enalaprilat")),
    has_other_ace       = any(compound %in% c("ramipril", "quinapril", "moexipril", "trandolapril")),
    has_despro          = any(compound %in% despro_compounds), .groups = "drop")

# Count files
a <- sum(ace_flags$has_enalapril_group &  ace_flags$has_despro)
b <- sum(ace_flags$has_enalapril_group & !ace_flags$has_despro)
c <- sum(ace_flags$has_other_ace       &  ace_flags$has_despro)
d <- sum(ace_flags$has_other_ace       & !ace_flags$has_despro)

# create a matrix based on counted files
cont_tab <- matrix(
  c(a, b, c, d), nrow = 2, byrow = TRUE,
  dimnames = list(
    DrugGroup = c("Enalapril(+Enalaprilat)", "Other_ACEI"),
    Despro    = c("Present", "Absent")))

cont_tab

# Fisher exact test
fisher_res <- fisher.test(cont_tab)

# odds ratio
fisher_res$estimate
# p value
fisher_res$p.value
# 95% CI
fisher_res$conf.int


# Build matrix

# restrict df to filename and compound
presence <- df %>%
  distinct(filename_level, compound) %>%
  dplyr::mutate(value = 1) %>%
  pivot_wider(names_from = compound, values_from = value, values_fill = 0)

# extract sample x compound matrix
mat <- presence %>% select(-filename_level)
compound_names <- colnames(mat)

cooccurrence <- matrix(0, nrow = length(compound_names), ncol = length(compound_names), dimnames = list(compound_names, compound_names))

for(i in seq_along(compound_names)) {
  for(j in seq_along(compound_names)) {
    cooccurrence[i, j] <- sum(mat[[i]] == 1 & mat[[j]] == 1)}}

cooccurrence <- as.data.frame(cooccurrence)
co_table <- as.data.frame(as.table(as.matrix(cooccurrence))) %>%
  dplyr::rename(compound_A = Var1, compound_B = Var2, cooccurrence = Freq) %>%
  rowwise() %>%
  dplyr::mutate(
    A_total = sum(mat[[compound_A]]),
    B_total = sum(mat[[compound_B]]),
    union_total = sum(mat[[compound_A]] | mat[[compound_B]]),
    pct_of_union = ifelse(union_total == 0, NA, cooccurrence / union_total * 100),
    pct_of_A = ifelse(A_total == 0, NA, cooccurrence / A_total * 100),
    pct_of_B = ifelse(B_total == 0, NA, cooccurrence / B_total * 100)) %>%
  ungroup()

co_table

# Create presence/absence table
presence <- df %>%
  distinct(filename_level, compound) %>%
  dplyr::mutate(value = 1) %>%
  pivot_wider(names_from = compound, values_from = value, values_fill = 0)

# Rename columns
presence <- presence %>%
  dplyr::rename(
    "desprolyl-enalaprilat" = enalaprilat_proline,
    "desprolyl-enalapril"   = enalapril_proline,
    "enalaprilat"           = Enalaprilat)

# Define molecule list with updated names
ace_molecules <- c("enalapril", "enalaprilat","ramipril", "quinapril", "trandolapril", "moexipril","desprolyl-enalapril", "desprolyl-enalaprilat")

# Select updated columns
plot_data <- presence %>% dplyr::select(all_of(ace_molecules))

# Upset plot
upset_plot <- ComplexUpset::upset(
  plot_data,
  ace_molecules,
  base_annotations = list('Intersection size' = intersection_size()),
  width_ratio = 0.3) +
  theme(
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.text.y  = element_text(color = "black"))

upset_plot
#ggsave("UpSet_plot_enalapril_drug.pdf", upset_plot, width = 8, height = 5, dpi = 900)
#getwd()
