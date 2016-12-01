# ==============================================================================
# run the biomarker TMLE across all of the genes
# ==============================================================================
library(foreach)
library(parallel)
library(doParallel)
library(ProjectTemplate)
load.project()
source("./lib/biomarkerTMLE.R")

# libraries for the g-fit and Q-fit steps of the TMLE
g.lib  <- c("SL.stepAIC")
Q.lib <- list("SL.glm", "SL.bayesglm")
#Q.lib <- list("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.nnet")


# ==============================================================================
# perform multi-level TMLE estimation (for all columns/genes)
# ==============================================================================
if (tail(strsplit(Sys.info()["nodename"], "[.]")$nodename, n = 1) == "edu") {
  doParallel::registerDoParallel(12)
} else {
  doParallel::registerDoParallel(4)
}

genomicATE <- foreach(gene = 1:ncol(Y2), .combine = cbind) %dopar% {
  print(paste("estimating ATE for", gene, "of", ncol(Y), ". Gene ID:",
              geneIDs[gene]))
  out <- biomarkerTMLE(Y = Y[, gene],
                       W = W,
                       A = A,
                       a = 1:length(unique(A)),
                       g.lib = g.lib,
                       Q.lib = Q.lib,
                       family = "gaussian"
                      )
}

biomarkerATE <- as.data.frame(t(genomicATE))
rownames(biomarkerATE) <- geneIDs
colnames(biomarkerATE) <- subjIDs
