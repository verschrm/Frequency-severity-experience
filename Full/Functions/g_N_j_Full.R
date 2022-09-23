##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the gradient vector of the
### expected complete log-likelihood of the number of 
### claims of profile j in the frequency-severity model 
### based on a Bayesian Hidden Markov Model.

### The function accepts the following arguments:

# theta_N_j    = Parameter vector for number of claims
#                of profile j;
# Prob_Z_hat_j = Estimated profile assignment 
#                probabilities for profile j;
# Counts       = Number of claims;
# Expo         = Exposure of policies in years;
# RF_N         = Risk factors for the claim counts;
# LB_eps       = Numerical precision for lower bounds;
# penalty      = Penalty for squared L2 norm of theta_N.

# Note that RF_N must contain numerical values and 
# that any categorical risk factors should thus have 
# already been one-hot encoded. Moreover, note that 
# theta_N_j assumes that the first ncol(RF_N) elements 
# represent delta_A^(j) and the last element a_U^(j).

##### Function g_N_j_Full() #####
g_N_j_Full <- function(theta_N_j, Prob_Z_hat_j, Counts, Expo, RF_N, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  delta_A_j <- theta_N_j[1:ncol(RF_N)]
  a_U_j <- max(tail(theta_N_j, 1), LB_eps)
  
  ### Scaled prior mean predictor
  Expo_lambda_j <- pmin(pmax(Expo * as.vector(exp(RF_N %*% delta_A_j)), LB_eps), 1e+15)
  
  ### Gradient components
  g_j_1 <- a_U_j * (Counts - Expo_lambda_j) / (a_U_j + Expo_lambda_j) * RF_N
  g_j_2 <- digamma(Counts + a_U_j) - digamma(a_U_j) + log(a_U_j) + 1 - 
    log(a_U_j + Expo_lambda_j) - (Counts + a_U_j) / (a_U_j + Expo_lambda_j)
  
  ### Gradient vector of log-likelihood
  g_j <- as.vector(Prob_Z_hat_j %*% cbind(g_j_1, g_j_2))
  g_j <- g_j - penalty * 2 * c(delta_A_j, a_U_j)
  
  ### Return output
  return(g_j)
}