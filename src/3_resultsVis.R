# visualization of results produced from statistical analyses

library(NMF)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("Rushmore", 100, type = "continuous")


# make heatmap of genes showing differential expression
top = 25
topgenes <- tt_diff_FDR %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice(1:top)
biomarkerATE <- biomarkerATE_diff %>%
  dplyr::filter(which(rownames(biomarkerATE_diff) %in% topgenes$geneID))

annot <- data.frame(Treatment = ifelse(design$Tx == 0, "Control", "Exposed"))
rownames(annot) <- colnames(biomarkerATE)
nmf.options(grid.patch = TRUE)

pdf(file = paste0(proj_dir, "/graphs/topGenesHeatmap.pdf"))
aheatmap(as.matrix(biomarkerATE),
         scale = "row",
         Rowv = TRUE,
         Colv = NULL,
         annCol = annot,
         annColors = "Set2",
         main = "Heatmap of Top 25 Genes (FDR-controlled)"
        )
dev.off()


# make histograms of raw and adjusted p-values and save to `graphs/` directory
pdf(file = paste0(proj_dir, "/graphs/limma_rawPval.pdf"))
ggplot(tt_diff, aes(P.Value)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.05) +
  ggtitle(paste("Histogram of raw p-values \n (applying Limma shrinkage to",
                "TMLE results)")) +
  xlab("magnitude of raw p-values") +
  scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) +
  theme_bw()
dev.off()

pdf(file = paste0(proj_dir, "/graphs/limma_adjPval.pdf"))
ggplot(tt_diff, aes(adj.P.Val)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.05) +
  ggtitle(paste("Histogram of FDR-corrected p-values (BH) \n (applying Limma",
                "shrinkage to TMLE results)")) +
  xlab("magnitude of BH-corrected p-values") +
  scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) +
  theme_bw()
dev.off()
