##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the expected complete
### log-likelihood of the claim sizes of profile j in 
### the frequency-severity model based on a Bayesian
### Hidden Markov Model.

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

##### Function Q_X_j_Full() #####
Q_X_j_Full <- function(theta_X_j, Prob_Z_hat_j, Sizes, RF_X, LB_eps = 1e-15,
                       penalty = 0) {
  ### Parameters
  phi_j <- max(theta_X_j[1], LB_eps)
  delta_B_j <- theta_X_j[2:(1 + ncol(RF_X))]
  a_V_j <- max(tail(theta_X_j, 1), 1 + LB_eps)
  
  ### Scaled prior mean predictor
  mu_j <- pmin(pmax(phi_j * as.vector(exp(RF_X %*% delta_B_j)), LB_eps), 1e+15)
  
  ### Natural logarithm of probabilities
  ln_Prob <- dbetapr(x = Sizes, shape1 = mu_j, shape2 = a_V_j,
                     scale = max((a_V_j - 1) / phi_j, LB_eps), log = TRUE)
  
  ### Log-likelihood
  Q_j <- as.numeric(Prob_Z_hat_j %*% ln_Prob)
  Q_j <- Q_j - penalty * (phi_j ^ 2 + sum(delta_B_j ^ 2) + a_V_j ^ 2)
  
  ### Return output
  return(Q_j)
}