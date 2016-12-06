# Targeted Minimum Loss-Based Estimation for the Analysis of Biomarkers
#
# This function performs influence curve-based estimation of the effect of an
# exposure on (gene) expression array values associated with a given biomarker,
# controlling for a user-specified set of baseline covariates
#
# Y (numeric vector) - a vector of array expression values for a single gene
# W (numeric matrix) - a matrix of covariates to be controlled in estimation
# A (numeric vector) - a discretized exposure vector whose effect on expression
#                      gene expression values is of interest
# a (numeric vector) - the levels of A against which comparisons are to be made
# g.lib (char vector) - library of learning algorithms to be used in ...
# Q.lib (char vector) - library of learning algorithms to be used in ...
# family (character) - specification of error family: "binomial" or "gaussian"
#
# NOTE - only compares the levels of A against the lowest level of A

biomarkerTMLE <- function(Y, W, A, a, param = "diff", family = "gaussian",
                          g.lib, Q.lib) {
  require(tmle)
  
  # check the case that Y is passed in as a column of a data.frame
  if (class(Y) == "data.frame") Y <- as.numeric(Y[, 1])
  if (class(A) == "data.frame") A <- as.numeric(A[, 1])
  if (length(a) == 1) warning("Comparisons should be made across levels of A")
  
  n_a = length(a)
  IC = NULL
  
  for(i in 1:n_a) {
    A_star = as.numeric(A == a[i])
    fit_tmle = tmle(Y = Y,
                    A = A_star,
                    W = W,
                    g.SL.library = g.lib,
                    Q.SL.library = Q.lib,
                    family = family,
                    verbose = FALSE
                   )
    g_0 = fit_tmle$g$g1W
    Qst_0 = fit_tmle$Qstar[, 2]
    EY_0 = mean(Qst_0)
    IC = cbind(IC, (A_star / g_0) * (Y - Qst_0) + Qst_0 - EY_0)
  }
  
  IC_diff = IC[, 2:n_a] - IC[, 1]
  
  if (param == "diff") {
    output <- IC_diff[, ncol(IC_diff)]
  } else {
    output <- IC[, n_a] # counterfactual difference between maximum Tx and not
  }
  return(output)
}
