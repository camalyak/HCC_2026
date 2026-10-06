#!/usr/bin/env Rscript


# Author: Kayla Mac
# Date: 10/6/26
# Data set: Plasmidsaurus HCC 
# Purpose: run DESeq2 on counts


# load necessary libraries
print("Loading libraries")

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

if (!require("DESeq2", quietly = TRUE)) 
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(DESeq2)

