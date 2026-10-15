#!/usr/bin/env Rscript


# Author: Kayla Mac
# Date: 10/6/26
# Data set: Plasmidsaurus HCC 
# Purpose: run DESeq2 on counts


# load necessary libraries

# load necessary libraries
print("Loading libraries")

suppressPackageStartupMessages({
  suppressWarnings({
    
    if (!require("DESeq2", quietly = TRUE)) {
      BiocManager::install("DESeq2", update = FALSE, ask = FALSE, quiet = TRUE)
    }
    
    if (!require("DEGreport", quietly = TRUE)) {
      BiocManager::install("DEGreport", update = FALSE, ask = FALSE, quiet = TRUE)
    }
    
    if (!require("EnhancedVolcano", quietly = TRUE)) {
      BiocManager::install("EnhancedVolcano", update = FALSE, ask = FALSE, quiet = TRUE)
    }
    
    if (!require("pacman", quietly = TRUE)) {
      install.packages("pacman", repos = "https://r-project.org", quiet = TRUE)
    }
    
    pacman::p_load(tidyverse, ggplot2, stringr, janitor, dplyr, DESeq2,
                   EnhancedVolcano)
    
  })
})


resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))

# creating volcano plots for diet, tumor, and aso affect on gene expression
res_diet_shr <- lfcShrink(dds, contrast = c("Diet", "HFD", "Chow"), type = "ashr")

lab <- mcols(dds)$gene_name[match(rownames(res_diet_shr), rownames(dds))]

EnhancedVolcano(res_diet_shr,
                lab = lab,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-6, 6),
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                max.overlaps = 5,
                labSize = 4,
                colConnectors = "grey40",
                arrowheads = FALSE)
