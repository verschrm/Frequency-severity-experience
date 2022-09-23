##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the Spearman's rho 
### and Kendall's tau implied by the Bayesian HMM.
### The function accepts the following arguments:

# K       = Number of latent risk profiles;
# size    = Number of successes for NB distribution;
# prob    = Probability of success for NB distribution;
# shape1  = First shape parameter of GB2 distribution;
# shape2  = Second shape parameter of GB2 distribution;
# shape3  = Third shape parameter of GB2 distribution;
# scale   = Scale parameter of GB2 distribution;
# Prob_Z  = Profile assignment probabilities.

##### Function rho_tau() #####
rho_tau <- function(K, size, prob, shape1, shape2, 
                    shape3, scale, Prob_Z) {
  ### Grid for u and v
  arg_grid <- seq(from = 1 / (2 * 100), to = 1 - 1 / (2 * 100), length.out = 100 + 1)
  
  ### Grid for number of claims
  N_grid <- seq(from = 0, to = 10, by = 1)
  
  ### Grid for claim size
  i <- 0
  X_grid <- list()
  for (v in arg_grid) { 
    i <- i + 1
    x_v <- as.vector(qgb2(prob = v, shape1 = shape1, scale = scale, shape2 = shape2, shape3 = shape3))
    X_grid[[i]] <- c(min(x_v), max(x_v))
  }
  X_grid <- seq(from = min(unlist(X_grid)), to = max(unlist(X_grid)), by = 1)
  
  ### Match each u to an n
  i <- 0
  cdf_N_grid <- list()
  for (n in N_grid) {
    i <- i + 1
    cdf_N_grid[[i]] <- sum(Prob_Z * sapply(1:K, function(j) 
      as.vector(pnbinom(q = n, size = size[j], prob = prob[j], lower.tail = TRUE, log.p = FALSE))))
  }
  cdf_N_grid <- matrix(unlist(cdf_N_grid), nrow = length(cdf_N_grid[[1]]), ncol = length(N_grid))
  i <- 0
  n_u <- list()
  for (u in arg_grid) {
    i <- i + 1
    n_u[[i]] <- N_grid[apply(abs(cdf_N_grid - u), 1, FUN = which.min)]
  }
  n_u <- matrix(unlist(n_u), nrow = length(n_u[[1]]), ncol = length(arg_grid))
  
  ### Match each v to an x
  i <- 0
  cdf_X_grid <- list()
  for (x in X_grid) {
    i <- i + 1
    cdf_X_grid[[i]] <- sum(Prob_Z * sapply(1:K, function(j) 
      as.vector(pgb2(x = x, shape1 = shape1, scale = scale[j], shape2 = shape2[j], shape3 = shape3[j]))))
  }
  cdf_X_grid <- matrix(unlist(cdf_X_grid), nrow = length(cdf_X_grid[[1]]), ncol = length(X_grid))
  i <- 0
  x_v <- list()
  for (v in arg_grid) {
    i <- i + 1
    x_v[[i]] <- X_grid[apply(abs(cdf_X_grid - v), 1, FUN = which.min)]
  }
  x_v <- matrix(unlist(x_v), nrow = length(x_v[[1]]), ncol = length(arg_grid))
  
  ### Midpoint rule for double integral
  Delta_arg <- arg_grid[2] - arg_grid[1]
  
  # CDF
  cdf_cop <- apply(Prob_Z * t(sapply(1:K, function(j)
    as.vector(pnbinom(q = n_u, size = size[j], prob = prob[j], lower.tail = TRUE, log.p = FALSE)) %o% as.vector(
      pgb2(x = x_v, shape1 = shape1, scale = scale[j], shape2 = shape2[j], shape3 = shape3[j])))), 2, sum)
  
  # PDF
  pdf_cop <- apply(Prob_Z * t(sapply(1:K, function(j)
    as.vector(dnbinom(x = n_u, size = size[j], prob = prob[j], log = FALSE)) %o% as.vector(
      dgb2(x = x_v, shape1 = shape1, scale = scale[j], shape2 = shape2[j], shape3 = shape3[j])))), 2, sum)
  pdf_cop <- pdf_cop / (apply(Prob_Z * t(sapply(1:K, function(j) as.vector(dnbinom(x = n_u, size = size[j], prob = prob[j], log = FALSE)))), 2, sum) %o% apply(
    Prob_Z * t(sapply(1:K, function(j) as.vector(dgb2(x = x_v, shape1 = shape1, scale = scale[j], shape2 = shape2[j], shape3 = shape3[j])))), 2, sum))
  pdf_cop[which(is.nan(pdf_cop))] <- 0
  pdf_cop_pmin <- pmin(pdf_cop, 1)
  
  # Dependence measure
  rho <- 12 * ((Delta_arg)^2 * sum(cdf_cop)) - 3
  tau <- 4 * ((Delta_arg)^2 * sum(cdf_cop * pdf_cop)) - 1
  
  return(c(rho, tau))
}