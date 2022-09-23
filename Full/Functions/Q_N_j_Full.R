##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the expected complete 
### log-likelihood of the number of claims of profile
### j in the frequency-severity model based on a 
### Bayesian Hidden Markov Model.

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

##### Function Q_N_j_Full() #####
Q_N_j_Full <- function(theta_N_j, Prob_Z_hat_j, Counts, Expo, RF_N, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  delta_A_j <- theta_N_j[1:ncol(RF_N)]
  a_U_j <- max(tail(theta_N_j, 1), LB_eps)
  
  ### Scaled prior mean predictor
  Expo_lambda_j <- pmin(pmax(Expo * as.vector(exp(RF_N %*% delta_A_j)), LB_eps), 1e+15)
  
  ### Natural logarithm of probabilities
  ln_Prob <- dnbinom(x = Counts, size = a_U_j, 
                     prob = pmax(a_U_j / (a_U_j + Expo_lambda_j), LB_eps), log = TRUE)
  
  ### Log-likelihood
  Q_j <- as.numeric(Prob_Z_hat_j %*% ln_Prob)
  Q_j <- Q_j - penalty * (sum(delta_A_j ^ 2) + a_U_j ^ 2)
  
  ### Return output
  return(Q_j)
}