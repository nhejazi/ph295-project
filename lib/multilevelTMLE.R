###########  Function below returns tmle estimates of risk differences of
###########  levels of A = a vs. the baseline (defined as smallest value of a)
tmle.multi = function(a = 1, A, W, Y, g.lib, Q.lib, family = "binomial") {
  library(tmle)
  # a is levels for which comparisons desired (assumes first level as baseline)
  IC = NULL; theta = NULL
  n_a = length(a)
  n_A = length(A)  # is there a reason that we define this?
  
  for(i in 1:n_a) {
    A_star = as.numeric(A == a[i])
    tmle.0 = tmle(Y = ytest,
                  A = A_star,
                  W = W,
                  g.SL.library = g.lib,
                  Q.SL.library = Q.lib,
                  family = family,
                  verbose = FALSE)
    g_0 = tmle.0$g$g1W
    Qst_0 = tmle.0$Qstar[, 2]
    Ey_0 = mean(Qst_0)
    IC = cbind(IC, (A_star / g_0) * (Y - Qst_0) + Qst_0 - Ey_0)
    theta = c(theta, mean(Qst_0))
  }
  
  ICd = IC[, 2:n_a] - IC[, 1] #last column for Limma
  psi = theta[2:n_a] - theta[1]
  psi = c(NA, psi)
  se = sqrt(apply(ICd, 2, var) / n_A)
  se = c(NA, se)
  lowr = psi - 1.96 * se
  uppr = psi + 1.96 * se
  CI95 = paste(round(lowr, 5), "-", round(uppr, 5), sep = "")
  pvalue = 2 * (1 - pnorm(abs(psi / se)))
  
  #alan_output = data.frame(a, theta, psi, se, CI95, pvalue)
  new_output = ICd[, 3]
  #output = list(new_output, alan_output)
  return(new_output)
}