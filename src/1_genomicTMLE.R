# ==============================================================================
# run the biomarker TMLE across all of the genes
# ==============================================================================
library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(ProjectTemplate)
load.project()
source("./lib/biomarkerTMLE.R")

# libraries for the g-fit and Q-fit steps of the TMLE
#g_lib_full  <- c("SL.glmnet", "SL.randomForest", "SL.nnet", "SL.polymars",
#                 "SL.mean", "SL.earth", "SL.gam", "SL.gbm")
#Q_lib_full <- c("SL.glmnet", "SL.randomForest", "SL.nnet", "SL.mean",
#                "SL.earth", "SL.gam")
g_lib_small  <- c("SL.glmnet", "SL.randomForest", "SL.nnet", "SL.polymars",
                  "SL.mean")
Q_lib_small <- c("SL.glmnet", "SL.randomForest", "SL.nnet", "SL.mean")

# ==============================================================================
# perform multi-level TMLE estimation (for all columns/genes)
# ==============================================================================
if (tail(strsplit(Sys.info()["nodename"], "[.]")$nodename, n = 1) == "edu") {
  doParallel::registerDoParallel(12)
} else {
  doParallel::registerDoParallel(4)
}

# computes the (rather complicated) parameter defined as the difference in the
# blips in Tx effects (difference between counterfactual max, min Tx effects)
genomicATE_diff <- foreach(gene = 1:ncol(Y_medNorm), .combine = cbind) %dopar% {
  print(paste("estimating ATE for", gene, "of", ncol(Y_medNorm), ". Gene ID:",
              geneIDs[gene]))
  out <- biomarkerTMLE(Y = Y_medNorm[, gene],
                       W = W,
                       A = A,
                       a = 1:length(unique(A)),
                       g.lib = g_lib_small,
                       Q.lib = Q_lib_small,
                       family = "gaussian"
                      )
}

biomarkerATE_diff <- as.data.frame(t(genomicATE_diff))
rownames(biomarkerATE_diff) <- geneIDs
colnames(biomarkerATE_diff) <- as.character(subjIDs)

data.table::fwrite(x = data.table(biomarkerATE_diff),
                   file = paste0(data_dir, "/IC_diff_libsmall.csv"))
