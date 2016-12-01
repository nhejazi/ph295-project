# this script takes the output generated from the multilevel (or regular) TMLE
# estimates for each gene, and implements the empirical Bayes / linear modeling
# method of Limma (basically, just call Limma on the IC-transformed estimates
# for each gene (that is, a gene x subject TMLE effects matrix))

library(limma)
library(ggplot2)
library(wesanderson)


# make the Limma design matrix and fit linear models to each gene
design <- as.data.frame(cbind(rep(1, nrow(Y)), as.numeric(A == min(unique(A)))))
colnames(design) <- c("intercept", "Tx")

fit <- lmFit(biomarkerATE, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, adjust.method = "BH", sort.by = "none",
               number = Inf)


# make histograms of raw and adjusted p-values and save to `graphs/` directory
pal <- wes_palette("Darjeeling", 100, type = "continuous")

pdf(paste0(proj_dir, "/graphs/limma_rawPval.pdf"))
ggplot(tt, aes(P.Value)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.01) +
  ggtitle("Histogram of raw p-values \n (Limma applied to TMLE)") +
  xlab("raw p-value magnitude") + scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) + xlim(0, 1) + theme_minimal()
dev.off()

pdf(paste0(proj_dir, "/graphs/limma_adjPval.pdf"))
ggplot(tt, aes(adj.P.Val)) +
  geom_histogram(aes(y = ..count.., fill = ..count..), colour = "white",
                 na.rm = TRUE, binwidth = 0.01) +
  ggtitle("Histogram of BH p-values (FDR) \n (Limma applied to TMLE)") +
  xlab("BH p-value magnitude") + scale_fill_gradientn("Count", colors = pal) +
  guides(fill = guide_legend(title = NULL)) + theme_minimal()
dev.off()