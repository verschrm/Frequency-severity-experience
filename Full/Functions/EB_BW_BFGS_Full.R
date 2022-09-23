##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function performs an empirical Bayes
### Baum-Welch algorithm with the BFGS method
### to maximize the frequency-severity log-likelihood 
### based on a Bayesian Hidden Markov Model. 
### It returns (i) the final parameter estimates 
### together with their sample variance estimates, 
### (ii) the posterior imputed profile assignment
### probabilities, (iii) the evolution of the 
### incomplete, or observed, log-likelihood and 
### (iii) the number of iterations in the algorithm.

### The function accepts the following arguments:

# theta_Z_init = Initial parameter matrix for theta_Z;
# theta_N_init = Initial parameter matrix for theta_N;
# theta_X_init = Initial parameter matrix for theta_X;
# K            = The number of latent risk profiles
#                to consider;
# Counts       = Number of claims;
# Expo         = Exposure of policies in years;
# Sizes        = Individual claim sizes;
# Ind_i_t      = Indicator of the policy underlying 
#                a claim;
# Ind_t_Ti     = Observation period and maximum 
#                number of periods observed;
# RF_N         = Risk factors for the claim counts;
# RF_X         = Risk factors for the claim sizes;
# Tol_abs_BW   = Absolute tolerance level of the
#                Baum-Welch algorithm;
# Tol_rel_BW   = Relative tolerance level of the
#                Baum-Welch algorithm;
# Iter_BW      = Maximum number of EM iterations;
# Tol_abs_M    = Absolute tolerance level of M-step;
# Tol_rel_M    = Relative tolerance level of M-step;
# Iter_M       = Maximum number of M-steps;
# Scaled       = Whether to scale the Baum-Welch 
#                algorithm (1) or not (0);
# LB_eps       = Numerical precision for lower bounds;
# penalty      = Penalty for squared L2 norm of 
#                theta_N and theta_X.

# Note that RF_N and RF_X must contain numerical
# values and that any categorical risk factors
# should thus have already been one-hot encoded.
# Moreover, note that theta_N assumes that the
# first rows represent delta_A and the last row
# a_U, and that theta_X assumes that the first
# row represent phi, the next rows delta_B and 
# the last row a_V. Finally, scaling of the 
# Baum-Welch algorithm usually increases its 
# numerical stability, especially in case of a 
# large number of periods.

##### Function EB_BW_BFGS_Full() #####
EB_BW_BFGS_Full <- function(theta_Z_init, theta_N_init, theta_X_init, 
                            K, Counts, Expo, Sizes, Ind_i_t, Ind_t_Ti,
                            RF_N, RF_X, Tol_abs_BW, Tol_rel_BW, Iter_BW, 
                            Tol_abs_M, Tol_rel_M, Iter_M, Scaled = 1, 
                            LB_eps = 1e-15, penalty = 0) {
  ### Ensure required packages are loaded
  Packages <- c('dplyr', 'maxLik', 'extraDistr', 'MASS')
  invisible(lapply(Packages, library, character.only = TRUE))
  
  ### Apply scaled version if invalid value is supplied
  Scaled <- ifelse(Scaled %in% c(0, 1), Scaled, 1)
  
  ### Initialization
  # Functions for computing the log-likelihood of the claim counts and sizes
  Q_N_j_Full <- dget('.../Full/Functions/Q_N_j_Full.R')
  Q_X_j_Full <- dget('.../Full/Functions/Q_X_j_Full.R')
  # Functions for computing the gradient vector of the log-likelihoods
  g_N_j_Full <- dget('.../Full/Functions/g_N_j_Full.R')
  g_X_j_Full <- dget('.../Full/Functions/g_X_j_Full.R')
  # Functions for computing the Hessian matrix of the log-likelihoods
  H_N_j_Full <- dget('.../Full/Functions/H_N_j_Full.R')
  H_X_j_Full <- dget('.../Full/Functions/H_X_j_Full.R')
  # Initialization of parameter values
  theta_Z <- theta_Z_init
  theta_N <- theta_N_init
  theta_X <- theta_X_init
  # Pre-allocate memory for variance estimates
  V_theta_Z <- NA * theta_Z
  V_theta_N <- NA * theta_N
  V_theta_X <- NA * theta_X
  # Tolerance levels
  Tol_abs <- 1
  Tol_rel <- 1
  # Number of iterations
  Iter <- 1
  Iter_N_r <- rep(0, K)
  Iter_X_r <- rep(0, K)
  # Incomplete, or observed, log-likelihood value
  ell_r <- 0
  # Convergence code
  Conv_N_r <- rep(0, K)
  Conv_X_r <- rep(0, K)
  
  ### Perform the Baum-Welch algorithm if the number of latent risk profiles is larger than 1
  if (K > 1) {
    ### Initial probabilities of observations
    # Parameters
    delta_A <- theta_N[1:ncol(RF_N), ]
    a_U <- as.vector(tail(theta_N, 1))
    phi <- as.vector(theta_X[1, ])
    delta_B <- theta_X[2:(1 + ncol(RF_X)), ]
    a_V <- as.vector(tail(theta_X, 1))
    # Transition probability matrix
    Prob_Z <- theta_Z
    # Determine scaled prior mean predictors
    Expo_lambda <- pmin(pmax(Expo * exp(RF_N %*% delta_A), LB_eps), 1e+15)
    mu <- pmin(pmax(t(phi * t(exp(RF_X %*% delta_B))), LB_eps), 1e+15)
    # Calculate claim count probability matrix
    Prob_N <- sapply(1:K, function(j) 
      dnbinom(x = Counts, size = a_U[j], 
              prob = pmax(a_U[j] / (a_U[j] + Expo_lambda[, j]), LB_eps), 
              log = FALSE))
    # Calculate claim size probability matrix
    Prob_X <- sapply(1:K, function(j)
      dbetapr(x = Sizes, shape1 = mu[, j], shape2 = a_V[j], 
              scale = pmax((a_V[j] - 1) / phi[j], LB_eps), log = FALSE))
    # Determine joint response, or emission, probabilities
    Prob_X_Total <- as.matrix(aggregate(Prob_X ~ Ind_i_t, FUN = prod))
    Prob_N_X <- Prob_N
    Prob_N_X[Prob_X_Total[, 1], ] <- Prob_N_X[Prob_X_Total[, 1], ] * Prob_X_Total[, -1]
    Prob_N_X <- pmax(Prob_N_X, .Machine$double.xmin)
    
    ### Initial probabilities of observation sequences
    # Compute forward filtering probabilities
    F_Prob <- matrix(1, nrow = nrow(RF_N), ncol = K)
    F_Prob[which(Ind_t_Ti[, 1] == 1), ] <- t(Prob_Z[1, ] * t(Prob_N_X[which(Ind_t_Ti[, 1] == 1), ]))
    # Scaled version for numerical stability
    if (Scaled == 1) {
      F_Scale <- as.vector(F_Prob[, 1])
      F_Scale[which(Ind_t_Ti[, 1] == 1)] <- apply(F_Prob[which(Ind_t_Ti[, 1] == 1), ], 1, sum)
      F_Prob[which(Ind_t_Ti[, 1] == 1), ] <- F_Prob[which(Ind_t_Ti[, 1] == 1), ] /
        F_Scale[which(Ind_t_Ti[, 1] == 1)]
    }
    for (t in 2:max(Ind_t_Ti[, 2])) {
     F_Prob[which(Ind_t_Ti[, 1] == t), ] <- Prob_N_X[which(Ind_t_Ti[, 1] == t), ] * 
        (F_Prob[which((Ind_t_Ti[, 1] == (t - 1)) & (Ind_t_Ti[, 2] >= t)), ] %*% Prob_Z[-1, ])
      # Scaled version for numerical stability
      if (Scaled == 1) {
        if (length(which(Ind_t_Ti[, 1] == t)) == 1) {
          F_Scale[which(Ind_t_Ti[, 1] == t)] <- sum(F_Prob[which(Ind_t_Ti[, 1] == t), ])
        } else {
          F_Scale[which(Ind_t_Ti[, 1] == t)] <- apply(F_Prob[which(Ind_t_Ti[, 1] == t), ], 1, sum)
        }
        F_Prob[which(Ind_t_Ti[, 1] == t), ] <- F_Prob[which(Ind_t_Ti[, 1] == t), ] /
          F_Scale[which(Ind_t_Ti[, 1] == t)]
      }
    }
    # Compute backward smoothing probabilities
    B_Prob <- matrix(1, nrow = nrow(RF_N), ncol = K)
    # Scaled version for numerical stability
    if (Scaled == 1) {
      B_Prob[which(Ind_t_Ti[, 1] == max(Ind_t_Ti[, 2])), ] <- 1 /
        F_Scale[which(Ind_t_Ti[, 1] == max(Ind_t_Ti[, 2]))]
    }
    for (t in (max(Ind_t_Ti[, 2]) - 1):1) {
      B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] <-
        (B_Prob[which(Ind_t_Ti[, 1] == (t + 1)), ] * Prob_N_X[which(Ind_t_Ti[, 1] == (t + 1)), ]) %*%
        t(Prob_Z[-1, ])
      # Scaled version for numerical stability
      if (Scaled == 1) {
        B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] <- 
          B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] / 
          F_Scale[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t))]
      }
    }
    
    ### Initial E-step for profile assignment probabilities
    # Impute posterior expected profile assignment probabilities
    Prob_Z_hat <- pmin(F_Prob * B_Prob, .Machine$double.xmax)
    Prob_Z_hat <- Prob_Z_hat / apply(Prob_Z_hat, 1, sum)
    
    ### Initial E-step for profile assignment probabilities starting from risk profile h
    # Impute posterior expected profile assignment probabilities starting from profile h
    Prob_Z_h_hat <- sapply(1:K, function(h) 
      F_Prob[which(Ind_t_Ti[, 1] < Ind_t_Ti[, 2]), h] * t(Prob_Z[1 + h, ] * t(
        B_Prob[which(Ind_t_Ti[, 1] > 1), ] * Prob_N_X[which(Ind_t_Ti[, 1] > 1), ])), simplify = 'array')
    Prob_Scale <- apply(Prob_Z_h_hat, 1, sum)
    Prob_Z_h_hat <- Prob_Z_h_hat / Prob_Scale
    
    ### Initial incomplete, or observed, log-likelihood
    if (Scaled == 0) {
      ell_r[1] <- sum(log(apply(F_Prob[which(Ind_t_Ti[, 1] == Ind_t_Ti[, 2]), ], 1, sum)))
    # Scaled version for numerical stability
    } else {
      ell_r[1] <- sum(log(F_Scale))
    }
    
    ### Print performance of initial solution (iteration 0)
    print(sprintf(paste('Iteration: %4i, Log-likelihood: %9.4f'), 0, ell_r[1]))
    
    ### Iterate until convergence or for Iter_BW iterations
    while ((Tol_abs > Tol_abs_BW) & (Tol_rel > Tol_rel_BW) & (Iter <= Iter_BW)) {
      ### M-step for initial profile assignment transitions
      theta_Z[1, ] <- pmax(apply(Prob_Z_hat[which(Ind_t_Ti[, 1] == 1), ], 2, mean), LB_eps)
      theta_Z[1, ] <- theta_Z[1, ] / sum(theta_Z[1, ])
      
      ### M-steps based on profile assignment probabilities
      # Pre-allocate memory for output
      Iter_N_r <- rbind(Iter_N_r, 0)
      Iter_X_r <- rbind(Iter_X_r, 0)
      Conv_N_r <- rbind(Conv_N_r, NA)
      Conv_X_r <- rbind(Conv_X_r, NA)
      for (j in 1:K) {
        ### M-step for claim counts conditional on risk profile j
        # Optimize complete expected log-likelihood using the BFGS method
        M_Q <- maxBFGS(fn = function(theta) Q_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat[, j], 
                                                       Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps,
                                                       penalty = penalty),
                       grad = function(theta) g_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat[, j], 
                                                         Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps,
                                                         penalty = penalty), 
                       hess = function(theta) H_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat[, j], 
                                                         Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps,
                                                         penalty = penalty), 
                       start = theta_N[, j], control = list(tol = Tol_abs_M, reltol = Tol_rel_M, 
                                                            iterlim = Iter_M, printLevel = 0), 
                       finalHessian = FALSE, constraints = list(ineqA = t(
                         c(rep(0, length(theta_N[, j]) - 1), 1)), ineqB = 0))
        # Store the resulting solution
        Iter_N_r[Iter + 1, j] <- as.numeric(M_Q$iterations[1])
        Conv_N_r[Iter + 1, j] <- M_Q$code
        theta_N[, j] <- M_Q$estimate
        # Ensure constraints are not violated
        theta_N[nrow(theta_N), j] <- pmax(theta_N[nrow(theta_N), j], LB_eps)
        
        ### M-step for individual claim sizes conditional on risk profile j
        # Optimize complete expected log-likelihood using the BFGS method
        M_Q <- maxBFGS(fn = function(theta) Q_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t, j], 
                                                       Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, 
                                                       penalty = penalty),
                       grad = function(theta) g_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t, j], 
                                                         Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, 
                                                         penalty = penalty), 
                       hess = function(theta) H_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t, j], 
                                                         Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, 
                                                         penalty = penalty), 
                       start = theta_X[, j], control = list(tol = Tol_abs_M, reltol = Tol_rel_M, 
                                                            iterlim = Iter_M, printLevel = 0), 
                       finalHessian = FALSE, constraints = list(ineqA = rbind(
                         c(1, rep(0, length(theta_X[, j]) - 1)), c(rep(0, length(theta_X[, j]) - 1), 1)),
                         ineqB = c(0, -1)))
        # Store the resulting solution
        Iter_X_r[Iter + 1, j] <- as.numeric(M_Q$iterations[1])
        Conv_X_r[Iter + 1, j] <- M_Q$code
        theta_X[, j] <- M_Q$estimate
        # Ensure constraints are not violated
        theta_X[1, j] <- pmax(theta_X[1, j], LB_eps)
        theta_X[nrow(theta_X), j] <- pmax(theta_X[nrow(theta_X), j], 1 + LB_eps)
      }
      
      ### M-steps based on profile assignment probabilities starting from profile h
      theta_Z[-1, ] <- pmax(t(apply(Prob_Z_h_hat, c(2, 3), sum)) / apply(Prob_Z_h_hat, 3, sum), LB_eps)
      theta_Z[-1, ] <- theta_Z[-1, ] / apply(theta_Z[-1, ], 1, sum)
      
      ### Update probabilities of observations
      # Parameters
      delta_A <- theta_N[1:ncol(RF_N), ]
      a_U <- as.vector(tail(theta_N, 1))
      phi <- as.vector(theta_X[1, ])
      delta_B <- theta_X[2:(1 + ncol(RF_X)), ]
      a_V <- as.vector(tail(theta_X, 1))
      # Transition probability matrix
      Prob_Z <- theta_Z
      # Determine scaled prior mean predictors
      Expo_lambda <- pmin(pmax(Expo * exp(RF_N %*% delta_A), LB_eps), 1e+15)
      mu <- pmin(pmax(t(phi * t(exp(RF_X %*% delta_B))), LB_eps), 1e+15)
      # Calculate claim count probability matrix
      Prob_N <- sapply(1:K, function(j) 
        dnbinom(x = Counts, size = a_U[j], 
                prob = pmax(a_U[j] / (a_U[j] + Expo_lambda[, j]), LB_eps), 
                log = FALSE))
      # Calculate claim size probability matrix
      Prob_X <- sapply(1:K, function(j)
        dbetapr(x = Sizes, shape1 = mu[, j], shape2 = a_V[j], 
                scale = pmax((a_V[j] - 1) / phi[j], LB_eps), log = FALSE))
      # Determine joint response, or emission, probabilities
      Prob_X_Total <- as.matrix(aggregate(Prob_X ~ Ind_i_t, FUN = prod))
      Prob_N_X <- Prob_N
      Prob_N_X[Prob_X_Total[, 1], ] <- Prob_N_X[Prob_X_Total[, 1], ] * Prob_X_Total[, -1]
      Prob_N_X <- pmax(Prob_N_X, .Machine$double.xmin)
      
      ### Update probabilities of observation sequences
      # Compute forward filtering probabilities
      F_Prob <- matrix(1, nrow = nrow(RF_N), ncol = K)
      F_Prob[which(Ind_t_Ti[, 1] == 1), ] <- t(Prob_Z[1, ] * t(Prob_N_X[which(Ind_t_Ti[, 1] == 1), ]))
      # Scaled version for numerical stability
      if (Scaled == 1) {
        F_Scale <- as.vector(F_Prob[, 1])
        F_Scale[which(Ind_t_Ti[, 1] == 1)] <- apply(F_Prob[which(Ind_t_Ti[, 1] == 1), ], 1, sum)
        F_Prob[which(Ind_t_Ti[, 1] == 1), ] <- F_Prob[which(Ind_t_Ti[, 1] == 1), ] /
          F_Scale[which(Ind_t_Ti[, 1] == 1)]
      }
      for (t in 2:max(Ind_t_Ti[, 2])) {
        F_Prob[which(Ind_t_Ti[, 1] == t), ] <- Prob_N_X[which(Ind_t_Ti[, 1] == t), ] * 
          (F_Prob[which((Ind_t_Ti[, 1] == (t - 1)) & (Ind_t_Ti[, 2] >= t)), ] %*% Prob_Z[-1, ])
        # Scaled version for numerical stability
        if (Scaled == 1) {
          if (length(which(Ind_t_Ti[, 1] == t)) == 1) {
            F_Scale[which(Ind_t_Ti[, 1] == t)] <- sum(F_Prob[which(Ind_t_Ti[, 1] == t), ])
          } else {
            F_Scale[which(Ind_t_Ti[, 1] == t)] <- apply(F_Prob[which(Ind_t_Ti[, 1] == t), ], 1, sum)
          }
          F_Prob[which(Ind_t_Ti[, 1] == t), ] <- F_Prob[which(Ind_t_Ti[, 1] == t), ] /
            F_Scale[which(Ind_t_Ti[, 1] == t)]
        }
      }
      # Compute backward smoothing probabilities
      B_Prob <- matrix(1, nrow = nrow(RF_N), ncol = K)
      # Scaled version for numerical stability
      if (Scaled == 1) {
        B_Prob[which(Ind_t_Ti[, 1] == max(Ind_t_Ti[, 2])), ] <- 1 /
          F_Scale[which(Ind_t_Ti[, 1] == max(Ind_t_Ti[, 2]))]
      }
      for (t in (max(Ind_t_Ti[, 2]) - 1):1) {
        B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] <-
          (B_Prob[which(Ind_t_Ti[, 1] == (t + 1)), ] * Prob_N_X[which(Ind_t_Ti[, 1] == (t + 1)), ]) %*%
          t(Prob_Z[-1, ])
        # Scaled version for numerical stability
        if (Scaled == 1) {
          B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] <- 
            B_Prob[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t)), ] / 
            F_Scale[which((Ind_t_Ti[, 1] == t) & (Ind_t_Ti[, 2] > t))]
        }
      }
      
      ### Update E-step for profile assignment probabilities
      # Impute posterior expected profile assignment probabilities
      Prob_Z_hat <- pmin(F_Prob * B_Prob, .Machine$double.xmax)
      Prob_Z_hat <- Prob_Z_hat / apply(Prob_Z_hat, 1, sum)
      
      ### Update E-step for profile assignment probabilities starting from risk profile h
      # Impute posterior expected profile assignment probabilities starting from profile h
      Prob_Z_h_hat <- sapply(1:K, function(h) 
        F_Prob[which(Ind_t_Ti[, 1] < Ind_t_Ti[, 2]), h] * t(Prob_Z[1 + h, ] * t(
          B_Prob[which(Ind_t_Ti[, 1] > 1), ] * Prob_N_X[which(Ind_t_Ti[, 1] > 1), ])), simplify = 'array')
      Prob_Scale <- apply(Prob_Z_h_hat, 1, sum)
      Prob_Z_h_hat <- Prob_Z_h_hat / Prob_Scale
      
      ### Update number of iterations, log-likelihood and tolerance levels
      Iter <- Iter + 1
      if (Scaled == 0) {
        ell_r <- c(ell_r, sum(log(apply(F_Prob[which(Ind_t_Ti[, 1] == Ind_t_Ti[, 2]), ], 1, sum))))
      # Scaled version for numerical stability
      } else {
        ell_r <- c(ell_r, sum(log(F_Scale)))
      }
      Tol_abs <- abs(ell_r[Iter] - ell_r[Iter - 1])
      Tol_rel <- abs(Tol_abs / ell_r[Iter - 1])
      
      ### Print performance of current iteration
      print(sprintf(paste('Iteration: %4i, Log-likelihood: %9.4f'), Iter - 1, ell_r[Iter]))
    }
    
    ### Final variance estimates
    # Initial profile assignment transitions
    V_theta_Z[1, ] <- theta_Z[1, ] / sum(Ind_t_Ti[, 1] == 1)
    # Risk profile j
    for (j in 1:K) {
      # Claim counts
      Hessian <- H_N_j_Full(theta_N_j = theta_N[, j], Prob_Z_hat_j = Prob_Z_hat[, j],
                            Counts = Counts, Expo = Expo, RF_N = RF_N)
      V_theta_N[, j] <- try(diag(solve(-Hessian)), silent = TRUE)
      if ((inherits(V_theta_N[, j], 'try-error')) | (sum(V_theta_N[, j] < 0) > 0)) {
        V_theta_N[, j] <- diag(ginv(-Hessian))
      }
      # Individual claim sizes
      Hessian <- H_X_j_Full(theta_X_j = theta_X[, j], Prob_Z_hat_j = Prob_Z_hat[Ind_i_t, j], 
                            Sizes = Sizes, RF_X = RF_X)
      V_theta_X[, j] <- try(diag(solve(-Hessian)), silent = TRUE)
      if ((inherits(V_theta_X[, j], 'try-error')) | (sum(V_theta_X[, j] < 0) > 0)) {
        V_theta_X[, j] <- diag(ginv(-Hessian))
      }
    }
    # Successive profile assignment probabilities starting from risk profile h
    V_theta_Z[-1, ] <- theta_Z[-1, ] / apply(Prob_Z_h_hat, 3, sum)
    
    ### Final posterior profile assignment probabilities for risk premium calculation
    Prob_Z_post_RP <- F_Prob / Prob_N_X
    Prob_Z_post_RP <- Prob_Z_post_RP / apply(Prob_Z_post_RP, 1, sum)
    
  ### Else perform the BFGS method directly on the entire portfolio
  } else {
    ### Artificial E-step
    Prob_Z_hat <- rep(1, nrow(RF_N))
    
    ### Initial log-likelihood
    ell_r[1] <- Q_N_j_Full(theta_N_j = theta_N, Prob_Z_hat_j = Prob_Z_hat, 
                           Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps, penalty = penalty) +
      Q_X_j_Full(theta_X_j = theta_X, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t], 
                 Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, penalty = penalty)
    
    ### Print performance of initial solution (iteration 0)
    print(sprintf(paste('Iteration: %4i, Log-likelihood: %9.4f'), 0, ell_r[1]))
    
    ### M-step for claim counts
    # Pre-allocate memory for output
    ell_r <- c(ell_r, 0)
    Iter_N_r <- c(Iter_N_r, 0)
    Conv_N_r <- c(Conv_N_r, NA)
    # Optimize log-likelihood using the BFGS method
    M_Q <- maxBFGS(fn = function(theta) Q_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat, 
                                                   Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps, 
                                                   penalty = penalty),
                   grad = function(theta) g_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat, 
                                                     Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps, 
                                                     penalty = penalty), 
                   hess = function(theta) H_N_j_Full(theta_N_j = theta, Prob_Z_hat_j = Prob_Z_hat, 
                                                     Counts = Counts, Expo = Expo, RF_N = RF_N, LB_eps = LB_eps, 
                                                     penalty = penalty), 
                   start = theta_N, control = list(tol = Tol_abs_M, reltol = Tol_rel_M, iterlim = Iter_M, 
                                                   printLevel = 0), finalHessian = TRUE, 
                   constraints = list(ineqA = t(c(rep(0, length(theta_N) - 1), 1)), ineqB = 0))
    # Store the resulting solution
    ell_r[2] <- M_Q$maximum
    Iter_N_r[2] <- as.numeric(M_Q$iterations[1])
    Conv_N_r[2] <- M_Q$code
    theta_N <- M_Q$estimate
    # Ensure constraints are not violated
    theta_N[nrow(theta_N)] <- max(theta_N[nrow(theta_N)], LB_eps)
    # Final variance estimates
    V_theta_N <- try(diag(solve(-M_Q$hessian)), silent = TRUE)
    if ((inherits(V_theta_N, 'try-error')) | (sum(V_theta_N < 0) > 0)) {
      V_theta_N <- diag(ginv(-M_Q$hessian))
    }
    
    ### M-step for individual claim sizes
    # Pre-allocate memory for output
    Iter_X_r <- c(Iter_X_r, 0)
    Conv_X_r <- c(Conv_X_r, NA)
    # Optimize log-likelihood using the BFGS method
    M_Q <- maxBFGS(fn = function(theta) Q_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t], 
                                                   Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, penalty = penalty),
                   grad = function(theta) g_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t], 
                                                     Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, penalty = penalty), 
                   hess = function(theta) H_X_j_Full(theta_X_j = theta, Prob_Z_hat_j = Prob_Z_hat[Ind_i_t], 
                                                     Sizes = Sizes, RF_X = RF_X, LB_eps = LB_eps, penalty = penalty), 
                   start = theta_X, control = list(tol = Tol_abs_M, reltol = Tol_rel_M, iterlim = Iter_M, 
                                                   printLevel = 0), finalHessian = TRUE, 
                   constraints = list(ineqA = rbind(c(1, rep(0, length(theta_X) - 1)), 
                                                    c(rep(0, length(theta_X) - 1), 1)), ineqB = c(0, -1)))
    # Store the resulting solution
    ell_r[2] <- ell_r[2] + M_Q$maximum
    Iter_X_r[2] <- as.numeric(M_Q$iterations[1])
    Conv_X_r[2] <- M_Q$code
    theta_X <- M_Q$estimate
    # Ensure constraints are not violated
    theta_X[1] <- max(theta_X[1], LB_eps)
    theta_X[nrow(theta_X)] <- max(theta_X[nrow(theta_X)], 1 + LB_eps)
    # Final variance estimates
    V_theta_X <- try(diag(solve(-M_Q$hessian)), silent = TRUE)
    if ((inherits(V_theta_X, 'try-error')) | (sum(V_theta_X < 0) > 0)) {
      V_theta_X <- diag(ginv(-M_Q$hessian))
    }
    
    ### Artificial posterior profile assignment probabilities for risk premium calculation
    Prob_Z_post_RP <- rep(1, nrow(RF_N))
    
    ### Update number of iterations
    Iter <- Iter + 1
    
    ### Print performance of current iteration
    print(sprintf(paste('Iteration: %4i, Log-likelihood: %9.4f'), 1, ell_r[2]))
  }
  
  ### Return output
  return(list('theta_Z' = theta_Z, 'V_theta_Z' = V_theta_Z,
              'theta_N' = theta_N, 'V_theta_N' = V_theta_N,
              'theta_X' = theta_X, 'V_theta_X' = V_theta_X,
              'Ind_Z_hat' = Prob_Z_hat, 'Prob_Z_post_RP' = Prob_Z_post_RP,
              'ell_BW' = ell_r, 'Iter_BW' = Iter,
              'Iter_M_N' = Iter_N_r, 'Conv_M_N' = Conv_N_r,
              'Iter_M_X' = Iter_X_r, 'Conv_M_X' = Conv_X_r))
}