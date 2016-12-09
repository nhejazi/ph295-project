# this script takes the output generated from the multilevel (or regular) TMLE
# estimates for each gene, and implements the empirical Bayes / linear modeling
# method of Limma (basically, just call Limma on the IC-transformed estimates
# for each gene (that is, a gene x subject TMLE effects matrix))

library(limma)

# make the Limma design matrix and fit linear models to each gene
design <- as.data.frame(cbind(rep(1, nrow(Y)), as.numeric(A == max(unique(A)))))
colnames(design) <- c("intercept", "Tx")

fit_diff <- lmFit(biomarkerATE_diff, design)
fit_diff <- eBayes(fit_diff)
tt_diff <- topTable(fit_diff, coef = 2, adjust.method = "BH",
                    sort.by = "none", number = Inf)
tt_diff$geneID <- rownames(biomarkerATE_diff)
data.table::fwrite(x = data.table(tt_diff),
                   file = paste0(proj_dir, "/results/topTableLimma.csv"))


# save table results of genes showing differential expression below FDR cutoff
FDRcutoff = 0.05
tt_diff_FDR <- tt_diff %>%
  subset(adj.P.Val < FDRcutoff)

data.table::fwrite(x = data.table(tt_diff),
                   file = paste0(proj_dir, "/results/topTableFDRLimma.csv"))
