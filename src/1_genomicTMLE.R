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
g_lib  <- c("SL.stepAIC", "SL.gbm", "SL.glmnet", "SL.loess", "SL.randomForest",
            "SL.nnet", "SL.knn", "SL.earth", "SL.polymars")
#Q_lib <- list("SL.glm", "SL.bayesglm")
#Q.lib <- list("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.nnet")


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
genomicATE_diff <- foreach(gene = 1:3, .combine = cbind) %dopar% {
  print(paste("estimating ATE for", gene, "of", ncol(Y_medNorm), ". Gene ID:",
              geneIDs[gene]))
  out <- biomarkerTMLE(Y = Y_medNorm[, gene],
                       W = W,
                       A = A,
                       a = 1:length(unique(A)),
                       g.lib = g_lib,
                       Q.lib = g_lib,
                       family = "gaussian"
                      )
}

biomarkerATE_diff <- as.data.frame(t(genomicATE_diff))
rownames(biomarkerATE_diff) <- geneIDs
colnames(biomarkerATE_diff) <- as.character(subjIDs)

data.table::fwrite(x = data.table(biomarkerATE_diff),
                   file = paste0(data_dir, "/ICestimates_diff.csv"))

# computes the parameter defined as the counterfactual difference in receiving
# the maximum Tx level versus not...
genomicATE_max <- foreach(gene = 1:ncol(Y), .combine = cbind) %dopar% {
  print(paste("estimating ATE for", gene, "of", ncol(Y_medNorm), ". Gene ID:",
              geneIDs[gene]))
  out <- biomarkerTMLE(Y = Y_medNorm[, gene],
                       W = W,
                       A = A,
                       a = 1:length(unique(A)),
                       g.lib = g.lib,
                       Q.lib = Q.lib,
                       family = "gaussian",
                       param = "max"
  )
}

biomarkerATE_max <- as.data.frame(t(genomicATE_max))
rownames(biomarkerATE_max) <- geneIDs
colnames(biomarkerATE_max) <- as.character(subjIDs)

data.table::fwrite(x = data.table(biomarkerATE_max),
                   file = paste0(data_dir, "/ICestimates_max.csv"))
