##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the gradient vector of the
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

##### Function g_X_j_Full() #####
g_X_j_Full <- function(theta_X_j, Prob_Z_hat_j, Sizes, RF_X, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  phi_j <- max(theta_X_j[1], LB_eps)
  delta_B_j <- theta_X_j[2:(1 + ncol(RF_X))]
  a_V_j <- max(tail(theta_X_j, 1), 1 + LB_eps)
  
  ### Scaled prior mean predictor
  mu_j <- pmin(pmax(phi_j * as.vector(exp(RF_X %*% delta_B_j)), LB_eps), 1e+15)
  
  ### Gradient components
  g_j_1 <- mu_j / phi_j * (digamma(mu_j + a_V_j) - digamma(mu_j) + log(phi_j) + 
                             1 + log(Sizes) - log(a_V_j - 1 + phi_j * Sizes)) -
    (mu_j + a_V_j) * Sizes / (a_V_j - 1 + phi_j * Sizes)
  g_j_2 <- mu_j * (digamma(mu_j + a_V_j) - digamma(mu_j) + log(phi_j) + log(Sizes) - 
                     log(a_V_j - 1 + phi_j * Sizes)) * RF_X
  g_j_3 <- digamma(mu_j + a_V_j) - digamma(a_V_j) + log(a_V_j - 1) + 
    1 / (a_V_j - 1) + 1 - log(a_V_j - 1 + phi_j * Sizes) - (mu_j + a_V_j) / 
    (a_V_j - 1 + phi_j * Sizes)
  
  ### Gradient vector of log-likelihood
  g_j <- as.vector(Prob_Z_hat_j %*% cbind(g_j_1, g_j_2, g_j_3))
  g_j <- g_j - penalty * 2 * c(phi_j, delta_B_j, a_V_j)
  
  ### Return output
  return(g_j)
}