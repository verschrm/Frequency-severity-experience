##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the gradient vector of the
### expected complete log-likelihood of the number of 
### claims in the frequency-severity model based on a 
### Bayesian Hidden Markov Model.

### The function accepts the following arguments:

# theta_N    = Parameter vector for number of claims;
# Prob_Z_hat = Estimated profile assignment 
#              probabilities for each profile;
# Counts     = Number of claims;
# Expo       = Exposure of policies in years;
# RF_N       = Risk factors for the claim counts;
# LB_eps     = Numerical precision for lower bounds;
# penalty    = Penalty for squared L2 norm of theta_N.

# Note that RF_N must contain numerical values and 
# that any categorical risk factors should thus have 
# already been one-hot encoded. Moreover, note that 
# theta_N assumes that its elements represent
# (delta_A, a_U, b_U).

##### Function g_N_Sparse() #####
g_N_Sparse <- function(theta_N, Prob_Z_hat, Counts, Expo, RF_N, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  K <- max(ncol(Prob_Z_hat), 1)
  delta_A <- as.vector(theta_N[1:ncol(RF_N)])
  a_U_mat <- matrix(pmax(theta_N[(ncol(RF_N) + 1):(ncol(RF_N) + K)], LB_eps),
                    nrow = nrow(RF_N), ncol = K, byrow = TRUE)
  b_U_mat <- matrix(pmax(theta_N[(ncol(RF_N) + 1 + K):(ncol(RF_N) + 2 * K)], LB_eps),
                    nrow = nrow(RF_N), ncol = K, byrow = TRUE)
  
  ### Scaled prior mean predictor
  Expo_lambda <- as.vector(pmin(pmax(Expo * exp(RF_N %*% delta_A), LB_eps), 1e+15))
  
  ### Gradient components
  g_1 <- as.vector(apply(apply(Prob_Z_hat * (Counts * b_U_mat - a_U_mat * Expo_lambda) / (
    b_U_mat + Expo_lambda), 1, sum) * RF_N, 2, sum))
  g_2 <- as.vector(apply(Prob_Z_hat * (digamma(Counts + a_U_mat) - digamma(a_U_mat) + 
                                         log(b_U_mat) - log(b_U_mat + Expo_lambda)), 2, sum))
  g_3 <- as.vector(apply(Prob_Z_hat * (a_U_mat * Expo_lambda - Counts * b_U_mat) / (
    b_U_mat * (b_U_mat + Expo_lambda)), 2, sum))
  
  ### Gradient vector of log-likelihood
  g <- as.vector(c(g_1, g_2, g_3))
  g <- g - penalty * 2 * c(delta_A, pmax(theta_N[(ncol(RF_N) + 1):(ncol(RF_N) + K)], LB_eps),
                           pmax(theta_N[(ncol(RF_N) + 1 + K):(ncol(RF_N) + 2 * K)], LB_eps))
  
  ### Return output
  return(g)
}