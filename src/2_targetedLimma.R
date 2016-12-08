# this script takes the output generated from the multilevel (or regular) TMLE
# estimates for each gene, and implements the empirical Bayes / linear modeling
# method of Limma (basically, just call Limma on the IC-transformed estimates
# for each gene (that is, a gene x subject TMLE effects matrix))

library(limma)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("Darjeeling", 100, type = "continuous")


# make the Limma design matrix and fit linear models to each gene
design <- as.data.frame(cbind(rep(1, nrow(Y)), as.numeric(A == max(unique(A)))))
colnames(design) <- c("intercept", "Tx")

fit_diff <- lmFit(biomarkerATE_diff, design)
fit_diff <- eBayes(fit_diff)
tt_diff <- topTable(fit_diff, coef = 2, adjust.method = "BH",
                    sort.by = "none", number = Inf)


# make histograms of raw and adjusted p-values and save to `graphs/` directory
pdf(paste0(proj_dir, "/graphs/limma_rawPval.pdf"))
ggplot(tt_diff, aes(P.Value)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.01) +
  ggtitle("Histogram of raw p-values \n (Limma applied to TMLE)") +
  xlab("raw p-value magnitude") + scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) + xlim(0, 1) + theme_minimal()
dev.off()

pdf(paste0(proj_dir, "/graphs/limma_adjPval.pdf"))
ggplot(tt_diff, aes(adj.P.Val)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.01) +
  ggtitle("Histogram of BH p-values (FDR) \n (Limma applied to TMLE)") +
  xlab("BH p-value magnitude") + scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) + xlim(0, 1) + theme_minimal()
dev.off()
