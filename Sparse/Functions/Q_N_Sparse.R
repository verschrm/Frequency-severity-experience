##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the expected complete 
### log-likelihood of the number of claims in the 
### frequency-severity model based on a Bayesian 
### Hidden Markov Model.

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

##### Function Q_N_Sparse() #####
Q_N_Sparse <- function(theta_N, Prob_Z_hat, Counts, Expo, RF_N, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  K <- max(ncol(Prob_Z_hat), 1)
  delta_A <- as.vector(theta_N[1:ncol(RF_N)])
  a_U <- as.vector(pmax(theta_N[(ncol(RF_N) + 1):(ncol(RF_N) + K)], LB_eps))
  b_U <- as.vector(pmax(theta_N[(ncol(RF_N) + 1 + K):(ncol(RF_N) + 2 * K)], LB_eps))
  
  ### Scaled prior mean predictor
  Expo_lambda <- as.vector(pmin(pmax(Expo * exp(RF_N %*% delta_A), LB_eps), 1e+15))
  
  ### Natural logarithm of probabilities
  ln_Prob <- sapply(1:K, function(j) 
    dnbinom(x = Counts, size = a_U[j], 
            prob = pmax(b_U[j] / (b_U[j] + Expo_lambda), LB_eps), log = TRUE))
  
  ### Log-likelihood
  Q <- as.numeric(sum(Prob_Z_hat * ln_Prob))
  Q <- Q - penalty * (delta_A ^ 2 + sum(a_U ^ 2) + sum(b_U ^ 2))
  
  ### Return output
  return(Q)
}