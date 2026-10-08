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

# filter low-count genes
print("Filtering low-count genes (at least one sample has 10 or more reads)")

keep <- rowSums(counts_clean >= 10) >= 1
print(paste("Genes before filtering:", nrow(counts_clean)))
print(paste("Genes after filtering: ", sum(keep)))

counts_clean <- counts_clean[keep, , drop = FALSE]


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

# running one-way ANOVA with multiple comparisons
print("running one-way ANOVA with multiple comparisons")

dds_grp <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData   = deseq_metadata,
  design    = ~ 0 + Group + Batch 
)

# fit the model
mcols(dds_grp) <- DataFrame(ensembl_ids)
dds_grp <- DESeq(dds_grp)
resultsNames(dds_grp)



raw    <- counts(dds_grp)
cd_all <- as.data.frame(colData(dds_grp))
cd_all$Mouse <- factor(cd_all$Mouse)

# one-way ANOVA + post-hoc, separately for each tissue
fit_tissue <- function(tis) {
  keep <- cd_all$Type == tis
  cd   <- cd_all[keep, ]
  cd$Group <- droplevels(factor(cd$Group))
  
  dds_t <- DESeqDataSetFromMatrix(raw[, keep], cd, design = ~ Group + Batch)
  mcols(dds_t) <- DataFrame(ensembl_ids)
  
  lrt  <- DESeq(dds_t, test = "LRT", reduced = ~ Batch)   # omnibus
  wald <- DESeq(dds_t)                                    # pairwise contrasts
  list(dds = wald, omni = results(lrt))
}

fits <- list(NonTumor = fit_tissue("NonTumor"), Tumor = fit_tissue("Tumor"))

# all 6 pairwise contrasts per tissue
res_within <- lapply(fits, function(f) {
  prs <- combn(levels(f$dds$Group), 2, simplify = FALSE)
  setNames(
    lapply(prs, function(p) results(f$dds, contrast = c("Group", p[2], p[1]))),
    sapply(prs, function(p) paste0(p[2], "_vs_", p[1]))
  )
})

# optional Bonferroni across all within-tissue contrasts (6 x 2 = 12)
n_contrasts <- sum(lengths(res_within))

# write CSVs
for (tis in names(fits)) {
  # omnibus
  write.csv(annotate_results(fits[[tis]]$omni, fits[[tis]]$dds),
            paste0("data/within_", tis, "_omnibus.csv"), row.names = FALSE)
  
  # pairwise
  for (nm in names(res_within[[tis]])) {
    r <- res_within[[tis]][[nm]]
    r$padj_all <- pmin(r$padj * n_contrasts, 1)
    res_within[[tis]][[nm]] <- r
    write.csv(annotate_results(r, fits[[tis]]$dds),
              paste0("data/within_", nm, ".csv"), row.names = FALSE)
  }
}

# how many genes are significant in each contrast?
sig_summary <- bind_rows(lapply(names(res_within), function(tis) {
  data.frame(
    tissue        = tis,
    contrast      = names(res_within[[tis]]),
    padj_0.05     = sapply(res_within[[tis]], function(r) sum(r$padj < 0.05, na.rm = TRUE)),
    padj_all_0.05 = sapply(res_within[[tis]], function(r) sum(r$padj_all < 0.05, na.rm = TRUE)),
    row.names = NULL
  )
}))
sig_summary
write.csv(sig_summary, "data/within_significance_summary.csv", row.names = FALSE)




# check the pairing first
table(cd_all$Mouse, cd_all$Type)

paired_res <- list()

for (d in c("Chow", "HFD")) for (a in c("Control", "ABAT")) {
  nm <- paste(d, a, sep = "_")
  in_cond <- cd_all$Diet == d & cd_all$ASO == a
  
  tab  <- table(droplevels(cd_all$Mouse[in_cond]), cd_all$Type[in_cond])
  both <- rownames(tab)[rowSums(tab > 0) == 2]     # mice with tumor AND non-tumor
  keep <- in_cond & cd_all$Mouse %in% both
  
  cd <- cd_all[keep, ]
  cd$Mouse <- droplevels(cd$Mouse)
  cd$Type  <- factor(cd$Type, levels = c("NonTumor", "Tumor"))
  
  sub <- DESeqDataSetFromMatrix(raw[, keep], cd, design = ~ Mouse + Type)
  mcols(sub) <- DataFrame(ensembl_ids)
  sub <- DESeq(sub)
  
  r <- results(sub, contrast = c("Type", "Tumor", "NonTumor"))
  paired_res[[nm]] <- r
  cat(nm, ":", length(both), "paired mice\n")
  
  write.csv(annotate_results(r, sub),
            paste0("data/paired_tumor_v_nontumor_", nm, ".csv"), row.names = FALSE)
}

# significant genes per condition
data.frame(
  condition = names(paired_res),
  padj_0.05 = sapply(paired_res, function(r) sum(r$padj < 0.05, na.rm = TRUE))
)
