# set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth")
library(tidyverse)
library(networkD3)
library(duckplyr)
library(htmlwidgets)

# import the FASST results from the python script
files <- list.files(
  path       = "/Users/vincentlamoureux/Desktop/antihypertensive_drug_FASST",
  pattern    = "matches\\.tsv$",
  full.names = TRUE)

# combine all files into a single dataframe
all_matches <- files |>
  purrr::map_dfr(~{
    compound_name <- basename(.x) |>
      str_remove("_matches\\.tsv$") |>
      str_remove("_")
    
    read_tsv(.x) |>
      dplyr::mutate(compound = compound_name)
  })

# create new columns and do some filtering
all_matches_datasets <- all_matches |> 
  dplyr::mutate(dataset = sub("^[^:]+:([^:]+):.*$", "\\1", USI)) |>
  dplyr::mutate(filename_level = str_remove(USI, ":scan:.*$")) |> 
  dplyr::filter(!compound == "enalaprilat_proline")

# remove non-human datasets
df <- all_matches_datasets |>
  dplyr::filter(!dataset %in% c("MSV000094528", "MSV000093469", "MSV000083593", "MSV000083772", "MSV000084908",
                         "MSV000085843", "MSV000087612", "MSV000087685", "MSV000092754", "MSV000094642",
                         "NORMAN-61a72f15-591b-5cf7-9a4e-7e26549f1564", "MSV000092487", "MSV000089190","MSV000096359",
                         "MSV000097575", "MSV000097976", "NORMAN-b3198128-7d7f-5d05-86a1-60da4b920140", "ST001715"))

#unique_datasets <- sort(unique(df$dataset))
#unique_datasets

# filter for the drug_metabolites
enalapril_proline <- df |>
  dplyr::filter(compound == "enalapril_proline") |>
  pull(filename_level) |>
  unique()

# compound that co-occur with enalapril_proline
co_occur <- df |>
  dplyr::filter(filename_level %in% enalapril_proline, compound != "enalapril_proline") |>
  pull(compound) |>
  unique()

# parent drug and metabolite should co-occur
parent_metabolite <- df |>
  dplyr::filter(compound != "enalapril_proline") |>
  distinct(dataset, filename_level, compound) |>
  count(dataset, compound, name = "value") |>
  dplyr::rename(source = dataset, target = compound)

# only co‑occurring compounds
only_co_occur <- df |>
  dplyr::filter(filename_level %in% enalapril_proline,
         compound %in% co_occur) |>
  distinct(filename_level, compound) |>
  count(compound, name = "value") |>
  dplyr::rename(source = compound) |>
  dplyr::mutate(target = "enalapril_proline")

# compounds that do NOT co‑occur with enalapril minus proline
compound_names   <- sort(unique(df$compound[df$compound != "enalapril_proline"]))
compounds_no_pro <- setdiff(compound_names, co_occur)

# gave a dummy node (no proline)
dummy <- tibble(
  source = compounds_no_pro,
  target = "no_proline",
  value  = 1)

# combine them
links <- bind_rows(parent_metabolite, only_co_occur, dummy)

# order the nodes
dataset_names <- sort(unique(df$dataset))
nodes_ordered <- c(dataset_names, compound_names, "enalapril_proline", "no_proline")
nodes <- data.frame(name = nodes_ordered, stringsAsFactors = FALSE)

# Links IDs
links <- links |>
  dplyr::mutate(
    IDsource = match(source, nodes$name) - 1,
    IDtarget = match(target, nodes$name) - 1)

# generate the plot
sn <- sankeyNetwork(
  Links       = links,
  Nodes       = nodes,
  Source      = "IDsource",
  Target      = "IDtarget",
  Value       = "value",
  NodeID      = "name",
  sinksRight  = TRUE,
  iterations  = 0,
  nodeWidth   = 40,
  nodePadding = 10,
  fontSize    = 12,
  width       = 500,
  height      = 900
)
sn
# export the plot to HTML
html_file <- "FigE_sankey.html"
saveWidget(sn, file = html_file, selfcontained = TRUE)
browseURL(html_file)
