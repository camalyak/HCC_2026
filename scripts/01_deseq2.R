#!/usr/bin/env Rscript


# Author: Kayla Mac
# Date: 10/6/26
# Data set: Plasmidsaurus HCC 
# Purpose: get Plasmidsaurus matrix and get rid of low yield samples


# load necessary libraries
print("Loading libraries")

suppressPackageStartupMessages({
  suppressWarnings({
    
    if (!require("org.Mm.eg.db", quietly = TRUE)) {
      BiocManager::install("org.Mm.eg.db", update = FALSE, ask = FALSE, quiet = TRUE)
    }
    if (!require("DESeq2", quietly = TRUE)) {
      BiocManager::install("DESeq2", update = FALSE, ask = FALSE, quiet = TRUE)
    }
    
    if (!require("pacman", quietly = TRUE)) {
      install.packages("pacman", repos = "https://r-project.org", quiet = TRUE)
    }
    
    pacman::p_load(tidyverse, ggplot2, stringr, janitor, dplyr)
    
  })
})


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
  mutate(Batch = case_when(
    nchar(as.character(Mouse)) == 3 ~ "A",
    nchar(as.character(Mouse)) == 4 ~ "B",
    TRUE ~ "Unknown"
  )) |> 
  mutate(Combined_ID = paste(Sample, Mouse, Type, Diet, ASO, sep = "_"))

counts_sample_ids <- as.numeric(sub("^[^_]+_([0-9]+)_.*$", "\\1", colnames(counts_clean)))

matched_indices <- match(counts_sample_ids, meta_combined$Sample)

colnames(counts_clean) <- meta_combined$Combined_ID[matched_indices]

meta_combined <- clean_names(meta_combined)
meta_combined <- meta_combined |>
  drop_na(sample, mouse)


# group by treatment group
print("Grouping by treatment group")

col_info <- data.frame(col_name = colnames(counts_clean), stringsAsFactors = FALSE) |> 
  separate(col_name, into = c("Sample", "Mouse", "Type", "Diet", "ASO"), sep = "_", remove = FALSE) |> 
  mutate(Batch = case_when(
    nchar(Mouse) == 3 ~ "A",
    nchar(Mouse) == 4 ~ "B",
    TRUE ~ "Unknown"
  ))

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

colnames(counts_reordered) <- paste0(col_info$Sample, "_", col_info$Mouse, "_", col_info$Group)


# add gene names that match the Ensembl ID
print("Matching Ensembl IDs to gene names")

original_ensembl <- rownames(counts_reordered)

clean_ensembl <- sub("\\..*$", "", original_ensembl)

gene_mappings <- mapIds(
  org.Mm.eg.db,
  keys = clean_ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

ensembl_ids <- data.frame(
  ensembl_id = original_ensembl,
  gene_name = unname(gene_mappings),
  stringsAsFactors = FALSE
) |> 
  mutate(gene_name = ifelse(is.na(gene_name), "Unknown", gene_name))

final_counts_with_genes <- cbind(ensembl_ids, as.data.frame(counts_reordered))

write.csv(final_counts_with_genes, "data/final_counts.csv", row.names = TRUE)
print("Dataframe successfully reordered and exported!")


# running DESeq2
deseq_metadata <- col_info |> 
  column_to_rownames("col_name")

deseq_counts <- as.matrix(counts_reordered)
deseq_counts <- round(deseq_counts)
storage.mode(deseq_counts) <- "integer"

print(paste("Metadata and Matrix match perfectly:", all(colnames(deseq_counts) == rownames(deseq_metadata))))
print(paste("Is matrix in integer format?:", is.integer(deseq_counts)))


print("running DESeq2: diet, aso, tumor effect on gene expression")
dds <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData = deseq_metadata,
  design = ~ Batch + Diet + Type + ASO
)

# Attach gene mapping tables inside the dds object metadata slot
mcols(dds) <- DataFrame(ensembl_ids)

dds <- DESeq(dds)

annotate_results <- function(res, dds) {
  anno <- as.data.frame(mcols(dds))[, c("ensembl_id", "gene_name")]
  stopifnot(identical(rownames(res), rownames(anno)))
  as.data.frame(res) |>
    cbind(anno) |>
    dplyr::relocate(ensembl_id, gene_name) |>
    dplyr::arrange(padj)
}

print("running diet effect on gene expression")
res_diet <- results(dds, contrast = c("Diet", "HFD", "Chow"))
summary(res_diet)
write.csv(annotate_results(res_diet, dds),
          "data/hfd_v_chow_annotated_results.csv", row.names = FALSE)

print("running aso effect on gene expression")
res_aso <- results(dds, contrast = c("ASO", "ABAT", "Control"))
summary(res_aso)
write.csv(annotate_results(res_aso, dds),
          "data/abat_v_ctrl_annotated_results.csv", row.names = FALSE)

print("running tumor effect on gene expression")
res_tum <- results(dds, contrast = c("Type", "Tumor", "NonTumor"))
summary(res_tum)
write.csv(annotate_results(res_tum, dds),
          "data/nt_v_t_annotated_results.csv", row.names = FALSE)







