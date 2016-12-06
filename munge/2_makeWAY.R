# create O = (W, A, Y) structure from cleaned up data


# W - age, sex, smoking
W <- data %>%
  dplyr::select(which(colnames(.) %in% c("age", "sex", "smoking")))


# A - benzene exposure (discretized)
A <- data %>%
  dplyr::select(which(colnames(.) %in% c("benzene")))
A <- A[, 1]


# Y - genes
Y <- data %>%
  dplyr::select(which(colnames(.) %ni% c("age", "sex", "smoking", "benzene",
                                         "id")))

geneIDs <- colnames(Y)

# median normalization
Y_medNorm <- as.data.frame(limma::normalizeBetweenArrays(Y, method = "scale"))

# (too simple of a) sanity check of whether Y includes array values
if(unique(lapply(Y_medNorm, class)) != "numeric") {
  print("Warning - values in Y do not appear to be gene expression measures...")
}
