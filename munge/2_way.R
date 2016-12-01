# create O = (W, A, Y) structure from cleaned up data

# W - age, sex, smoking
W <- data %>%
  select(which(colnames(.) %in% c("age", "sex", "smoking")))

# A - benzene exposure (discretized)
A <- data %>%
  select(which(colnames(.) %in% c("benzene")))

# Y - genes
Y <- data %>%
  select(which(colnames(.) %ni% c("age", "sex", "smoking", "benzene", "id")))

# (too simple of a) sanity check of whether Y includes array values
if(unique(lapply(Y, class)) != "numeric") {
  print("Warning - values in Y do not appear to be gene expression measures...")
}