##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the Hessian matrix of the
### expected complete log-likelihood of the claim sizes 
### in the frequency-severity model based on a Bayesian
### Hidden Markov Model.

### The function accepts the following arguments:

# theta_X    = Parameter vector for claim sizes;
# Prob_Z_hat = Estimated profile assignment 
#              probabilities for each profile;
# Sizes      = Individual claim sizes;
# RF_X       = Risk factors for the claim sizes;
# LB_eps     = Numerical precision for lower bounds;
# penalty    = Penalty for squared L2 norm of theta_X.

# Note that RF_X must contain numerical values and 
# that any categorical risk factors should thus have 
# already been one-hot encoded. Moreover, note that 
# theta_X assumes that its elements represent
# (phi, delta_B, a_V, b_V).

##### Function H_X_Sparse() #####
H_X_Sparse <- function(theta_X, Prob_Z_hat, Sizes, RF_X, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  K <- max(ncol(Prob_Z_hat), 1)
  phi <- as.numeric(max(theta_X[1], LB_eps))
  delta_B <- as.vector(theta_X[2:(1 + ncol(RF_X))])
  a_V_mat <- matrix(pmax(theta_X[(2 + ncol(RF_X)):(1 + ncol(RF_X) + K)], 1 + LB_eps),
                    nrow = nrow(RF_X), ncol = K, byrow = TRUE)
  b_V_mat <- matrix(pmax(theta_X[(2 + ncol(RF_X) + K):(1 + ncol(RF_X) + 2 * K)], LB_eps),
                    nrow = nrow(RF_X), ncol = K, byrow = TRUE)
  
  ### Scaled prior mean predictor
  mu <- as.vector(pmin(pmax(phi * exp(RF_X %*% delta_B), LB_eps), 1e+15))
  
  ### Hessian components
  H_1_1 <- as.numeric(sum(Prob_Z_hat / phi * (mu / phi * (mu * (
    trigamma(mu + a_V_mat) - trigamma(mu)) + 1 - phi * Sizes / (b_V_mat + phi * Sizes)) - 
      (b_V_mat * mu - a_V_mat * phi * Sizes) * Sizes / (b_V_mat + phi * Sizes) ^ 2)))
  H_1_2 <- as.vector(apply(mu / phi * apply(Prob_Z_hat * (mu * (
    trigamma(mu + a_V_mat) - trigamma(mu)) + digamma(mu + a_V_mat) - digamma(mu) + log(phi) + 
      1 + log(Sizes) - log(b_V_mat + phi * Sizes) - phi * Sizes / 
      (b_V_mat + phi * Sizes)), 1, sum) * RF_X, 2, sum))
  H_1_3 <- as.vector(apply(Prob_Z_hat * (mu / phi * trigamma(mu + a_V_mat) - Sizes / 
                                           (b_V_mat + phi * Sizes)), 2, sum))
  H_1_4 <- as.vector(apply(Prob_Z_hat / phi * (a_V_mat * phi * Sizes - b_V_mat * mu) /
                             (b_V_mat + phi * Sizes) ^ 2, 2, sum))
  H_2_2 <- t(RF_X) %*% (mu * apply(Prob_Z_hat * (mu * (
    trigamma(mu + a_V_mat) - trigamma(mu)) + digamma(mu + a_V_mat) - digamma(mu) + log(phi) + 
      log(Sizes) - log(b_V_mat + phi * Sizes)), 1, sum) * RF_X)
  H_2_3 <- t(RF_X) %*% (Prob_Z_hat * mu * trigamma(mu + a_V_mat))
  H_2_4 <- - t(RF_X) %*% (Prob_Z_hat * mu / (b_V_mat + phi * Sizes))
  H_3_3 <- diag(as.vector(apply(Prob_Z_hat * (trigamma(mu + a_V_mat) - trigamma(a_V_mat)), 
                                2, sum)), nrow = K)
  H_3_4 <- diag(as.vector(apply(Prob_Z_hat * phi * Sizes / (b_V_mat * (b_V_mat + phi * Sizes)), 
                                2, sum)), nrow = K)
  H_4_4 <- diag(as.vector(- apply(Prob_Z_hat * (a_V_mat * phi * Sizes * (
    phi * Sizes + 2 * b_V_mat) - mu * b_V_mat ^ 2) / (b_V_mat ^ 2 * (b_V_mat + phi * Sizes) ^ 2), 
    2, sum)), nrow = K)
  
  ### Hessian matrix of log-likelihood
  H <- rbind(c(H_1_1, H_1_2, H_1_3, H_1_4), 
             cbind(H_1_2, H_2_2, H_2_3, H_2_4), 
             cbind(H_1_3, t(H_2_3), H_3_3, H_3_4),
             cbind(H_1_4, t(H_2_4), t(H_3_4), H_4_4))
  H <- H - diag(2 * penalty, nrow = nrow(H), ncol = ncol(H))
  
  ### Return output
  return(H)
}