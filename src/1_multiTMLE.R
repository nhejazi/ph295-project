# ==============================================================================
# run multi-level TMLE
# ==============================================================================
# load helper TMLE function and use ProjectTemplate
library(ProjectTemplate)
load.project()
source("./lib/multilevelTMLE.R")

# set parameters
n.sample <- length(Y)
aa <- 1:4
g.lib  <- c("SL.stepAIC")
#Q.lib <- list("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.step.interaction",
#              "SL.nnet")
Q.lib <- list("SL.glm", "SL.bayesglm")


# ==============================================================================
# perform multi-level TMLE estimation (for one column only)
# THIS IS FOR TESTING WHETHER THE FUNCTION WORKS ONLY
# ==============================================================================
#out <- tmle.multi(aa,
#                  A = A[, 1], W = W, Y = Y$Y,
#                  g.lib, Q.lib, fam = "binomial")


# ==============================================================================
# perform multi-level TMLE estimation (for all columns)
# ==============================================================================
library(foreach); library(doParallel)
n_Cores <- detectCores()
registerDoParallel(n_Cores)

ATE_4vs1 <- foreach(gene = 1:ncol(A), .combine = rbind) %dopar% {
  print(paste("estimating casual effect", gene, "of", ncol(A), ". Gene ID:",
              colnames(A)[gene]))
  out <- tmle.multi(aa,
                    A = A[, gene], W = W, Y = Y$Y,
                    g.lib, Q.lib, fam = "binomial")
}
ATE_4vs1 <- as.data.frame(ATE_4vs1)