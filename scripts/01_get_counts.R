#!/usr/bin/env Rscript


# Author: Kayla Mac
# Date: 10/6/26
# Data set: Plasmidsaurus HCC 
# Purpose: get Plasmidsaurus matrix and get rid of low yield samples


# load necessary libraries
print("Loading libraries")

BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(tidyverse,
               ggplot2,
               stringr,
               janitor,
               dplyr,
               org.Hs.eg.db,
               AnnotationDbi)


# load raw matrix
print("Getting raw expression matrix")

mt <- read.table("plasmidsaurus_results/F3HP8S-expression-matrix.tsv",
                     header = TRUE,
                     row.names = 1,
                     fill = TRUE)


# get counts matrix
print("Getting counts matrix")

counts <- as.matrix(mt[, grepl("count", colnames(mt)), drop = FALSE])


# reorder matrix ascending
print("Sorting by ascending")

sorted_sample_names <- str_sort(colnames(counts), numeric = TRUE)
counts <- counts[, sorted_sample_names]


# remove unwanted samples
print("Getting rid of unwanted samples (low yield or wrong cancer marker)")

smpl_ls <- read.csv("sample_list/hcc_unwanted_sample_list.csv")
smpl_ls <- clean_names(smpl_ls)
smpl_ls <- remove_empty(smpl_ls, "cols")

smpl_ls <- smpl_ls |>
  drop_na(sample, mouse)

unwanted_samples <- smpl_ls$sample

extracted_numbers <- as.numeric(sub("^[^_]+_([0-9]+)_.*$", "\\1", colnames(counts)))
counts_clean <- counts[, !extracted_numbers %in% unwanted_samples]


# rename samples to treatment groups
print("Renaming samples to treatment groups")

meta <- read.csv("sample_list/hcc_wanted_sample_list.csv")

meta_combined <- meta |> 
  mutate(Combined_ID = paste(Sample, Mouse, Type, Diet, ASO, sep = "_"))

counts_sample_ids <- as.numeric(sub("^[^_]+_([0-9]+)_.*$", "\\1", colnames(counts_clean)))

matched_indices <- match(counts_sample_ids, meta_combined$Sample)

colnames(counts_clean) <- meta_combined$Combined_ID[matched_indices]


# group by treatment group
print("Grouping by treatment group")

col_info <- data.frame(col_name = colnames(counts_clean), stringsAsFactors = FALSE) |> 
  separate(col_name, into = c("Sample", "Mouse", "Type", "Diet", "ASO"), sep = "_", remove = FALSE)

desired_order <- c(
  "NonTumor_Chow_Control",
  "NonTumor_Chow_ABAT",
  "NonTumor_HFD_Control",
  "NonTumor_HFD_ABAT",
  "Tumor_Chow_Control",
  "Tumor_Chow_ABAT",
  "Tumor_HFD_Control",
  "Tumor_HFD_ABAT"
)


col_info <- col_info |> 
  mutate(Group = paste(Type, Diet, ASO, sep = "_")) |> 
  mutate(Group = factor(Group, levels = desired_order)) |> 
  arrange(Group, as.numeric(Sample))

counts_reordered <- counts_clean[, col_info$col_name]

colnames(counts_reordered) <- paste0("Sample_", col_info$Sample, "_", col_info$Group)


# add gene names that match the Ensembl ID
ensembl_ids <- data.frame(
  ensembl_id = original_ensembl,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)



















write.csv(counts_reordered, "data/reordered_treatment_counts.csv", row.names = TRUE)
print("Dataframe successfully reordered and exported!")






