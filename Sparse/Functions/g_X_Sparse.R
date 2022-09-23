##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the gradient vector of the
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

##### Function g_X_Sparse() #####
g_X_Sparse <- function(theta_X, Prob_Z_hat, Sizes, RF_X, LB_eps = 1e-15,
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
  
  ### Gradient components
  g_1 <- as.numeric(sum(Prob_Z_hat * (mu / phi * (
    digamma(mu + a_V_mat) - digamma(mu) + log(phi) + 1 + log(Sizes) -
      log(b_V_mat + phi * Sizes)) - (mu + a_V_mat) * Sizes / (b_V_mat + phi * Sizes))))
  g_2 <- as.vector(apply(mu * apply(Prob_Z_hat * (
    digamma(mu + a_V_mat) - digamma(mu) + log(phi) + log(Sizes) - 
      log(b_V_mat + phi * Sizes)), 1, sum) * RF_X, 2, sum))
  g_3 <- as.vector(apply(Prob_Z_hat * (
    digamma(mu + a_V_mat) - digamma(a_V_mat) + log(b_V_mat) - 
      log(b_V_mat + phi * Sizes)), 2, sum))
  g_4 <- as.vector(apply(Prob_Z_hat / b_V_mat * (a_V_mat * phi * Sizes - mu * b_V_mat) /
                           (b_V_mat + phi * Sizes), 2, sum))
  
  ### Gradient vector of log-likelihood
  g <- as.vector(c(g_1, g_2, g_3, g_4))
  g <- g - penalty * 2 * c(phi, delta_B, pmax(theta_X[(2 + ncol(RF_X)):(1 + ncol(RF_X) + K)], 1 + LB_eps),
                           pmax(theta_X[(2 + ncol(RF_X) + K):(1 + ncol(RF_X) + 2 * K)], LB_eps))
  
  ### Return output
  return(g)
}