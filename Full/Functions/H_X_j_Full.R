##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the Hessian matrix of the
### expected complete log-likelihood of the claim sizes 
### of profile j in the frequency-severity model based on 
### a Bayesian Hidden Markov Model.

### The function accepts the following arguments:

# theta_X_j    = Parameter vector for claim sizes
#                of profile j;
# Prob_Z_hat_j = Estimated profile assignment 
#                probabilities for profile j;
# Sizes        = Individual claim sizes;
# RF_X         = Risk factors for the claim sizes;
# LB_eps       = Numerical precision for lower bounds;
# penalty      = Penalty for squared L2 norm of theta_X.

# Note that RF_X must contain numerical values and 
# that any categorical risk factors should thus have 
# already been one-hot encoded. Moreover, note that 
# theta_X_j assumes that the first element represents 
# phi^(j), the next ncol(RF_X) elements delta_B^(j) and 
# the last element a_V^(j).

##### Function H_X_j_Full() #####
H_X_j_Full <- function(theta_X_j, Prob_Z_hat_j, Sizes, RF_X, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  phi_j <- max(theta_X_j[1], LB_eps)
  delta_B_j <- theta_X_j[2:(1 + ncol(RF_X))]
  a_V_j <- max(tail(theta_X_j, 1), 1 + LB_eps)
  
  ### Scaled prior mean predictor
  mu_j <- pmin(pmax(phi_j * as.vector(exp(RF_X %*% delta_B_j)), LB_eps), 1e+15)
  
  ### Hessian components
  H_j_1_1 <- mu_j / phi_j * (mu_j / phi_j * (
    trigamma(mu_j + a_V_j) - trigamma(mu_j)) + 1 / phi_j - 2 * Sizes / 
      (a_V_j - 1 + phi_j * Sizes)) + (mu_j + a_V_j) * 
    (Sizes / (a_V_j - 1 + phi_j * Sizes)) ^ 2
  H_j_1_2 <- mu_j / phi_j * (digamma(mu_j + a_V_j) - digamma(mu_j) + mu_j * (
    trigamma(mu_j + a_V_j) - trigamma(mu_j)) + log(phi_j) + 1 + log(Sizes) - 
      log(a_V_j - 1 + phi_j * Sizes) - phi_j * Sizes / 
      (a_V_j - 1 + phi_j * Sizes)) * RF_X
  H_j_1_3 <- (mu_j * trigamma(mu_j + a_V_j) - (
    mu_j * (a_V_j - 1) + phi_j * Sizes * (phi_j * Sizes - 1)) / 
      (a_V_j - 1 + phi_j * Sizes) ^ 2) / phi_j
  H_j_2_2 <- t(RF_X) %*% (Prob_Z_hat_j * mu_j * (
    digamma(mu_j + a_V_j) - digamma(mu_j) + mu_j * (
      trigamma(mu_j + a_V_j) - trigamma(mu_j)) + log(phi_j) + log(Sizes) - 
      log(a_V_j - 1 + phi_j * Sizes)) * RF_X)
  H_j_2_3 <- mu_j * (trigamma(mu_j + a_V_j) - 1 / 
                       (a_V_j - 1 + phi_j * Sizes)) * RF_X
  H_j_3_3 <- trigamma(mu_j + a_V_j) - trigamma(a_V_j) + (a_V_j - 2) / 
    (a_V_j - 1) ^ 2 - (a_V_j - 2 + 2 * phi_j * Sizes - mu_j) / 
    (a_V_j - 1 + phi_j * Sizes) ^ 2
  
  ### Hessian matrix of log-likelihood
  H_j_1 <- as.vector(Prob_Z_hat_j %*% cbind(H_j_1_1, H_j_1_2, H_j_1_3))
  H_j_3m1 <- as.vector(Prob_Z_hat_j %*% cbind(H_j_2_3, H_j_3_3))
  H_j <- cbind(H_j_1, 
               rbind(H_j_1[2:(1 + ncol(RF_X))], H_j_2_2, H_j_3m1[1:ncol(RF_X)]), 
               as.vector(c(tail(H_j_1, 1), H_j_3m1)))
  H_j <- H_j - diag(2 * penalty, nrow = nrow(H_j), ncol = ncol(H_j))
  
  ### Return output
  return(H_j)
}