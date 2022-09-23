##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function calculates the expected complete
### log-likelihood of the claim sizes in the 
### frequency-severity model based on a Bayesian
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

##### Function Q_X_Sparse() #####
Q_X_Sparse <- function(theta_X, Prob_Z_hat, Sizes, RF_X, LB_eps = 1e-15, 
                       penalty = 0) {
  ### Parameters
  K <- max(ncol(Prob_Z_hat), 1)
  phi <- as.numeric(max(theta_X[1], LB_eps))
  delta_B <- as.vector(theta_X[2:(1 + ncol(RF_X))])
  a_V <- as.vector(pmax(theta_X[(2 + ncol(RF_X)):(1 + ncol(RF_X) + K)], 1 + LB_eps))
  b_V <- as.vector(pmax(theta_X[(2 + ncol(RF_X) + K):(1 + ncol(RF_X) + 2 * K)], LB_eps))
  
  ### Scaled prior mean predictor
  alpha <- as.vector(pmin(pmax(phi * exp(RF_X %*% delta_B), LB_eps), 1e+15))
  
  ### Natural logarithm of probabilities
  ln_Prob <- sapply(1:K, function(j)
    dbetapr(x = Sizes, shape1 = alpha, shape2 = a_V[j], scale = max(b_V[j] / phi, LB_eps), 
            log = TRUE))
  
  ### Log-likelihood
  Q <- as.numeric(sum(Prob_Z_hat * ln_Prob))
  Q <- Q - penalty * (phi ^ 2 + sum(delta_B ^ 2) + sum(a_V ^ 2) + sum(b_V ^ 2))
  
  ### Return output
  return(Q)
}