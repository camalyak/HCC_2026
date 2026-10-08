#!/usr/bin/env Rscript


# Author: Kayla Mac
# Date: 10/8/26
# Data set: Plasmidsaurus HCC 
# Purpose: determine if batch and by aso differ 

vsd <- vst(dds_grp, blind = TRUE)

plotPCA(vsd, intgroup = c("Group"))
pca <- plotPCA(vsd, intgroup = c("Batch", "Type", "Diet", "ASO"), returnData = TRUE)
pca$name <- rownames(pca)
pv <- round(100 * attr(pca, "percentVar"))

p <- ggplot(pca, aes(PC1, PC2, color = Type, shape = Batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) +
  theme_bw()

ggsave("figures/pca/pca_type_batch_group.png", plot = p, width = 7, height = 5, dpi = 300)

p2 <- ggplot(pca, aes(PC1, PC2, color = ASO, shape = Batch)) +
  geom_point(size = 3) +
  facet_wrap(~ Type) +
  theme_bw()

ggsave("figures/pca/pca_type_batch.png", plot = p2, width = 7, height = 5, dpi = 300)

pcs <- prcomp(t(assay(vsd)[order(-rowVars(assay(vsd)))[1:500], ]))
plot(pcs$x[, 3], pcs$x[, 4], col = factor(vsd$Type), pch = as.numeric(factor(vsd$Batch)) + 15)