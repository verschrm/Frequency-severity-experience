##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Clear workspace #####
cat("\014")
rm(list = ls())
while (dev.cur()>1) dev.off()

#### Libraries #####
Packages <- c('dplyr', 'maxLik', 'extraDistr', 
              'MASS', 'gamlss', 'matrixStats')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load results from MMTPL_Preprocessing.R

##### Input settings #####
### Load function for Emprical Bayes Baum-Welch algorithm with BFGS method
EB_BW_BFGS_Full_Full <- dget('.../Full/Functions/EB_BaumWelch_BFGS_Full.R')

### Tolerance and iteration settings for Baum-Welch algorithm
Tol_abs_BW <- 1e-8
Tol_rel_BW <- 1e-8
Iter_BW <- 250
penalty <- 4e-6
# Increased precision
Tol_abs_BW_incr <- 1e-12
Tol_rel_BW_incr <- 1e-12
Iter_BW_incr <- 50000

### Tolerance and iteration setting for M-step with BFGS method
LB_eps <- 1e-15
Tol_abs_M <- 1e-8
Tol_rel_M <- 1e-8
Iter_M <- 100000
# Increased precision
Tol_abs_M_incr <- 1e-12
Tol_rel_M_incr <- 1e-12

##### Optimize (K=1) #####
### Initial parameter estimates for a single risk profile (K=1)
# Number of latent risk profiles
K <- 1
# Latent profile assignments
theta_Z_init_1 <- NA
NB_init <- glm.nb(formula = MTPL_N_Freq ~ -1 + MTPL_RF_Freq + offset(log(MTPL_Exp_Freq)), 
                  link = log, control = glm.control(epsilon = Tol_rel_M, maxit = Iter_M))
theta_N_init_1 <- c(as.vector(NB_init$coefficients), NB_init$theta)
# Individual claim sizes
G_init <- glm(formula = MTPL_X_Sev ~ -1 + MTPL_RF_Sev, family = Gamma(link = 'log'))
G_init$shape <- gamma.shape(G_init)
GB2_init <- gamlss(formula = MTPL_X_Sev ~ 1, sigma.formula = MTPL_X_Sev ~ 1,
                   nu.formula = MTPL_X_Sev ~ -1 + MTPL_RF_Sev, tau.formula = MTPL_X_Sev ~ 1,
                   family = GB2(mu.link = 'log', sigma.link = 'log',
                                nu.link = 'log', tau.link = 'log'),
                   mu.start = (2 - 1) / G_init$shape$alpha, sigma.start = 1,
                   nu.start = G_init$shape$alpha * as.numeric(G_init$fitted.values),
                   tau.start = 2, sigma.fix = TRUE,
                   control = gamlss.control(c.crit = Tol_abs_M, n.cyc = Iter_M, trace = FALSE), 
                   i.control = glim.control(cc = Tol_abs_M, cyc = Iter_M, bf.cyc = Iter_M, 
                                            bf.tol = Tol_abs_M))
theta_X_init_1 <- as.numeric((exp(GB2_init$tau.coefficients) - 1) / exp(GB2_init$mu.coefficients))
theta_X_init_1 <- c(theta_X_init_1, as.numeric(GB2_init$nu.coefficients[1] - log(theta_X_init_1)),
                    as.vector(GB2_init$nu.coefficients[-1]), as.numeric(exp(GB2_init$tau.coefficients)))

### Optimize for K=1 using the Empirical Bayes Baum-Welch algorithm with the BFGS method
Optim_1 <- EB_BW_BFGS_Full_Full(theta_Z_init = theta_Z_init_1, theta_N_init = theta_N_init_1,
                           theta_X_init = theta_X_init_1, K = K, Counts = MTPL_N_Freq, 
                           Expo = MTPL_Exp_Freq, Sizes = MTPL_X_Sev, Ind_i_t = MTPL_Ind_i_t, 
                           Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, RF_X = MTPL_RF_Sev, 
                           Tol_abs_BW = Tol_abs_BW_incr, Tol_rel_BW = Tol_rel_BW_incr, 
                           Iter_BW = Iter_BW_incr, Tol_abs_M = Tol_abs_M_incr, 
                           Tol_rel_M = Tol_rel_M_incr, Iter_M = Iter_M, Scaled = 1, 
                           LB_eps = LB_eps, penalty = penalty)

##### Optimize (K=2) #####
### Search grid for initial parameter values for two risk profiles (K=2)
# Number of latent risk profiles
K <- 2
K_distort_2 <- seq(from = qnorm(0.025), to = qnorm(0.975), length.out = K)
Dim_init_2 <- 2 * (11 - 1)
## Latent profile assignments
theta_Z_init_2 <- array(NA, dim = c(1 + K, K, Dim_init_2))
theta_Z_init_2[, , 1] <- t(sapply(1:(K + 1), function(j) rep(1 / K, K)))
theta_Z_init_2[, , 2] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 3] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 4] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 5] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 6] <- theta_Z_init_2[, , 1]
# -
theta_Z_init_2[, , 7] <- rbind(c(0.75, 0.25), c(0.75, 0.25), c(0.25, 0.75))
theta_Z_init_2[, , 8] <- theta_Z_init_2[, , 7]
theta_Z_init_2[, , 9] <- theta_Z_init_2[, , 7]
theta_Z_init_2[, , 10] <- theta_Z_init_2[, , 7]
theta_Z_init_2[, , 11] <- theta_Z_init_2[, , 7]
theta_Z_init_2[, , 12] <- theta_Z_init_2[, , 7]
theta_Z_init_2[, , 13] <- theta_Z_init_2[, , 7]
# -
theta_Z_init_2[, , 14] <- rbind(c(0.85, 0.15), c(0.95, 0.05), c(0.05, 0.95))
theta_Z_init_2[, , 15] <- theta_Z_init_2[, , 14]
theta_Z_init_2[, , 16] <- theta_Z_init_2[, , 14]
theta_Z_init_2[, , 17] <- theta_Z_init_2[, , 14]
theta_Z_init_2[, , 18] <- theta_Z_init_2[, , 14]
theta_Z_init_2[, , 19] <- theta_Z_init_2[, , 14]
theta_Z_init_2[, , 20] <- theta_Z_init_2[, , 14]
## Number of claims
theta_N_init_2 <- array(NA, dim = c(length(Optim_1$theta_N), K, Dim_init_2))
theta_N_init_2[, , 1] <- sapply(1:K, function(j) Optim_1$theta_N + K_distort_2[j] * 
                                  sqrt(pmax(Optim_1$V_theta_N, Tol_abs_M ^ 2)))
theta_N_init_2[length(Optim_1$theta_N), , 1] <- 
  theta_N_init_2[length(Optim_1$theta_N), c(2, 1), 1]
theta_N_init_2[, , 2] <- theta_N_init_2[, , 1]
theta_N_init_2[, , 3] <- cbind(c(2 ^ 2 * Optim_1$theta_N[1], 
                                 Optim_1$theta_N[-c(1, length(Optim_1$theta_N))], 
                                 2 * Optim_1$theta_N[length(Optim_1$theta_N)]),
                               c(0.5 ^ 2 * Optim_1$theta_N[1], 
                                 Optim_1$theta_N[-c(1, length(Optim_1$theta_N))],
                                 0.5 * Optim_1$theta_N[length(Optim_1$theta_N)]))
theta_N_init_2[, , 4] <- theta_N_init_2[, , 3]
theta_N_init_2[, , 5] <- Optim_1$theta_N * 
  cbind(ifelse(Optim_1$theta_N < 0, 2, 0.5), ifelse(Optim_1$theta_N < 0, 0.5, 2))
theta_N_init_2[, , 6] <- theta_N_init_2[, , 5]
# -
theta_N_init_2[, , 7] <- cbind(Optim_1$theta_N, Optim_1$theta_N)
theta_N_init_2[, , 8] <- theta_N_init_2[, , 1]
theta_N_init_2[, , 9] <- theta_N_init_2[, , 2]
theta_N_init_2[, , 10] <- theta_N_init_2[, , 3]
theta_N_init_2[, , 11] <- theta_N_init_2[, , 4]
theta_N_init_2[, , 12] <- theta_N_init_2[, , 5]
theta_N_init_2[, , 13] <- theta_N_init_2[, , 6]
# -
theta_N_init_2[, , 14] <- theta_N_init_2[, , 7]
theta_N_init_2[, , 15] <- theta_N_init_2[, , 1]
theta_N_init_2[, , 16] <- theta_N_init_2[, , 2]
theta_N_init_2[, , 17] <- theta_N_init_2[, , 3]
theta_N_init_2[, , 18] <- theta_N_init_2[, , 4]
theta_N_init_2[, , 19] <- theta_N_init_2[, , 5]
theta_N_init_2[, , 20] <- theta_N_init_2[, , 6]
# Individual claim sizes
theta_X_init_2 <- array(NA, dim = c(length(Optim_1$theta_X), K, Dim_init_2))
theta_X_init_2[, , 1] <- sapply(1:K, function(j) Optim_1$theta_X + K_distort_2[j] * 
                                  sqrt(pmax(Optim_1$V_theta_X, Tol_abs_M ^ 2)))
theta_X_init_2[c(1, length(Optim_1$theta_X)), , 1] <- 
  theta_X_init_2[c(1, length(Optim_1$theta_X)), c(2, 1), 1]
theta_X_init_2[, , 2] <- theta_X_init_2[, c(2, 1), 1]
theta_X_init_2[, , 3] <- cbind(c(c(2, 0.5) ^ 2 * Optim_1$theta_X[c(1, 2)], 
                                 Optim_1$theta_X[-c(1, 2, length(Optim_1$theta_X))],
                                 2 * Optim_1$theta_X[length(Optim_1$theta_X)]),
                               c(c(0.5, 2) ^ 2 * Optim_1$theta_X[c(1, 2)], 
                                 Optim_1$theta_X[-c(1, 2, length(Optim_1$theta_X))],
                                 0.5 * Optim_1$theta_X[length(Optim_1$theta_X)]))
theta_X_init_2[, , 4] <- theta_X_init_2[, c(2, 1), 3]
theta_X_init_2[, , 5] <- Optim_1$theta_X * 
  cbind(ifelse(Optim_1$theta_X < 0, 2, 0.5), ifelse(Optim_1$theta_X < 0, 0.5, 2))
theta_X_init_2[, , 6] <- theta_X_init_2[, c(2, 1), 5]
# -
theta_X_init_2[, , 7] <- cbind(Optim_1$theta_X, Optim_1$theta_X)
theta_X_init_2[, , 8] <- theta_X_init_2[, , 1]
theta_X_init_2[, , 9] <- theta_X_init_2[, , 2]
theta_X_init_2[, , 10] <- theta_X_init_2[, , 3]
theta_X_init_2[, , 11] <- theta_X_init_2[, , 4]
theta_X_init_2[, , 12] <- theta_X_init_2[, , 5]
theta_X_init_2[, , 13] <- theta_X_init_2[, , 6]
# -
theta_X_init_2[, , 14] <- theta_X_init_2[, , 7]
theta_X_init_2[, , 15] <- theta_X_init_2[, , 1]
theta_X_init_2[, , 16] <- theta_X_init_2[, , 2]
theta_X_init_2[, , 17] <- theta_X_init_2[, , 3]
theta_X_init_2[, , 18] <- theta_X_init_2[, , 4]
theta_X_init_2[, , 19] <- theta_X_init_2[, , 5]
theta_X_init_2[, , 20] <- theta_X_init_2[, , 6]
# Ensure constraints are met at initial values
theta_X_init_2[1, , ] <- pmax(theta_X_init_2[1, , ], 0.10)
theta_X_init_2[dim(theta_X_init_2)[1], , ] <- 
  pmax(theta_X_init_2[dim(theta_X_init_2)[1], , ], 1.10)
theta_N_init_2[dim(theta_N_init_2)[1], , ] <- 
  pmax(theta_N_init_2[dim(theta_N_init_2)[1], , ], 0.10)

### Try different starting values
ell_init_2 <- rep(NA, Dim_init_2)
for (init in 1:Dim_init_2) {
  ell_init_2[init] <- tail(EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_2[, , init], 
                                           theta_N_init = theta_N_init_2[, , init],
                                           theta_X_init = theta_X_init_2[, , init], K = K, 
                                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, 
                                           Sizes = MTPL_X_Sev, Ind_i_t = MTPL_Ind_i_t, 
                                           Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW, 
                                           Tol_rel_BW = Tol_rel_BW, Iter_BW = Iter_BW, 
                                           Tol_abs_M = Tol_abs_M, Tol_rel_M = Tol_rel_M, 
                                           Iter_M = Iter_M, Scaled = 1, LB_eps = LB_eps, 
                                           penalty = penalty)$ell_BW, 1)
  print(sprintf(paste('Initial run: %2i / %2i, Log-likelihood: %9.4f'), init, Dim_init_2, 
                ell_init_2[init]))
}
best_init_2 <- which.max(ell_init_2)

### Optimize for K=2 using the Empirical Bayes Baum-Welch algorithm with the BFGS method
Optim_2 <- EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_2[, , best_init_2], 
                           theta_N_init = theta_N_init_2[, , best_init_2],
                           theta_X_init = theta_X_init_2[, , best_init_2], K = K, 
                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, Sizes = MTPL_X_Sev, 
                           Ind_i_t = MTPL_Ind_i_t, Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW_incr, 
                           Tol_rel_BW = Tol_rel_BW_incr, Iter_BW = Iter_BW_incr, 
                           Tol_abs_M = Tol_abs_M_incr, Tol_rel_M = Tol_rel_M_incr, Iter_M = Iter_M, 
                           Scaled = 1, LB_eps = LB_eps, penalty = penalty)

##### Optimize (K=3) #####
### Search grid for initial parameter values for three risk profiles (K=3)
# Number of latent risk profiles
K <- 3
K_distort_3 <- seq(from = qnorm(0.025), to = qnorm(0.975), length.out = K)
Dim_init_3 <- 6 * (2 * 3 + 2)
## Latent profile assignments
theta_Z_init_3 <- array(NA, dim = c(1 + K, K, Dim_init_3))
theta_Z_init_3[, , 1] <- t(sapply(1:(K + 1), function(j) rep(1 / K, K)))
theta_Z_init_3[, , 2] <- theta_Z_init_3[, , 1]
theta_Z_init_3[, , 3] <- theta_Z_init_3[, , 2]
theta_Z_init_3[, , 4] <- theta_Z_init_3[, , 3]
theta_Z_init_3[, , 5] <- theta_Z_init_3[, , 4]
theta_Z_init_3[, , 6] <- theta_Z_init_3[, , 5]
theta_Z_init_3[, , 7] <- theta_Z_init_3[, , 6]
theta_Z_init_3[, , 8] <- theta_Z_init_3[, , 7]
theta_Z_init_3[, , 9] <- theta_Z_init_3[, , 8]
theta_Z_init_3[, , 10] <- theta_Z_init_3[, , 9]
theta_Z_init_3[, , 11] <- theta_Z_init_3[, , 10]
theta_Z_init_3[, , 12] <- theta_Z_init_3[, , 11]
# -
theta_Z_init_3[, , 13] <- rbind(c(0.50, 0.30, 0.20), c(0.50, 0.30, 0.20), 
                                c(0.30, 0.50, 0.20), c(0.20, 0.30, 0.50))
theta_Z_init_3[, , 14] <- theta_Z_init_3[, , 13]
theta_Z_init_3[, , 15] <- theta_Z_init_3[, , 14]
theta_Z_init_3[, , 16] <- theta_Z_init_3[, , 15]
theta_Z_init_3[, , 17] <- theta_Z_init_3[, , 16]
theta_Z_init_3[, , 18] <- theta_Z_init_3[, , 17]
theta_Z_init_3[, , 19] <- theta_Z_init_3[, , 18]
theta_Z_init_3[, , 20] <- theta_Z_init_3[, , 19]
theta_Z_init_3[, , 21] <- theta_Z_init_3[, , 20]
theta_Z_init_3[, , 22] <- theta_Z_init_3[, , 21]
theta_Z_init_3[, , 23] <- theta_Z_init_3[, , 22]
theta_Z_init_3[, , 24] <- theta_Z_init_3[, , 23]
theta_Z_init_3[, , 25] <- theta_Z_init_3[, , 24]
theta_Z_init_3[, , 26] <- theta_Z_init_3[, , 25]
theta_Z_init_3[, , 27] <- theta_Z_init_3[, , 26]
theta_Z_init_3[, , 28] <- theta_Z_init_3[, , 27]
theta_Z_init_3[, , 29] <- theta_Z_init_3[, , 28]
theta_Z_init_3[, , 30] <- theta_Z_init_3[, , 29]
# -
theta_Z_init_3[, , 31] <- rbind(c(0.75, 0.15, 0.10), c(0.80, 0.125, 0.075), 
                                c(0.125, 0.80, 0.075), c(0.075, 0.125, 0.80))
theta_Z_init_3[, , 32] <- theta_Z_init_3[, , 31]
theta_Z_init_3[, , 33] <- theta_Z_init_3[, , 32]
theta_Z_init_3[, , 34] <- theta_Z_init_3[, , 33]
theta_Z_init_3[, , 35] <- theta_Z_init_3[, , 34]
theta_Z_init_3[, , 36] <- theta_Z_init_3[, , 35]
theta_Z_init_3[, , 37] <- theta_Z_init_3[, , 36]
theta_Z_init_3[, , 38] <- theta_Z_init_3[, , 37]
theta_Z_init_3[, , 39] <- theta_Z_init_3[, , 38]
theta_Z_init_3[, , 40] <- theta_Z_init_3[, , 39]
theta_Z_init_3[, , 41] <- theta_Z_init_3[, , 40]
theta_Z_init_3[, , 42] <- theta_Z_init_3[, , 41]
theta_Z_init_3[, , 43] <- theta_Z_init_3[, , 42]
theta_Z_init_3[, , 44] <- theta_Z_init_3[, , 43]
theta_Z_init_3[, , 45] <- theta_Z_init_3[, , 44]
theta_Z_init_3[, , 46] <- theta_Z_init_3[, , 45]
theta_Z_init_3[, , 47] <- theta_Z_init_3[, , 46]
theta_Z_init_3[, , 48] <- theta_Z_init_3[, , 47]
## Number of claims
theta_N_init_3 <- array(NA, dim = c(nrow(Optim_2$theta_N), K, Dim_init_3))
theta_N_init_3[, , 1] <- cbind(c(2 ^ 2 * Optim_2$theta_N[1, 1], Optim_2$theta_N[
  -c(1, nrow(Optim_2$theta_N)), 1], 0.5 * Optim_2$theta_N[nrow(Optim_2$theta_N), 1]), NA, 
  c(0.5 ^ 2 * Optim_2$theta_N[1, 2], Optim_2$theta_N[-c(1, nrow(Optim_2$theta_N)), 2],
    2 * Optim_2$theta_N[nrow(Optim_2$theta_N), 2]))
theta_N_init_3[, 2, 1] <- apply(theta_N_init_3[, c(1, 3), 1], 1, mean)
theta_N_init_3[, , 2] <- theta_N_init_3[, , 1]
theta_N_init_3[, , 3] <- theta_N_init_3[, , 2]
theta_N_init_3[, , 4] <- theta_N_init_3[, , 3]
theta_N_init_3[, , 5] <- theta_N_init_3[, , 4]
theta_N_init_3[, , 6] <- theta_N_init_3[, , 5]
theta_N_init_3[, , 7] <- cbind(Optim_2$theta_N[, 1] * c(ifelse(Optim_2$theta_N[
  -nrow(Optim_2$theta_N), 1] < 0, 2, 0.5), 0.5), NA, Optim_2$theta_N[, 2] * c(ifelse(
    Optim_2$theta_N[-nrow(Optim_2$theta_N), 2] < 0, 0.5, 2), 2))
theta_N_init_3[, 2, 7] <- apply(theta_N_init_3[, c(1, 3), 7], 1, mean)
theta_N_init_3[, , 8] <- theta_N_init_3[, , 7]
theta_N_init_3[, , 9] <- theta_N_init_3[, , 8]
theta_N_init_3[, , 10] <- theta_N_init_3[, ,9]
theta_N_init_3[, , 11] <- theta_N_init_3[, , 10]
theta_N_init_3[, , 12] <- theta_N_init_3[, , 11]
# -
theta_N_init_3[, , 13] <- cbind(Optim_2$theta_N[, 1], apply(Optim_2$theta_N, 1, mean),
                                Optim_2$theta_N[, 2])
theta_N_init_3[, , 14] <- theta_N_init_3[, , 13]
theta_N_init_3[, , 15] <- theta_N_init_3[, , 14]
theta_N_init_3[, , 16] <- theta_N_init_3[, , 15]
theta_N_init_3[, , 17] <- theta_N_init_3[, , 16]
theta_N_init_3[, , 18] <- theta_N_init_3[, , 17]
theta_N_init_3[, , 19] <- theta_N_init_3[, , 1]
theta_N_init_3[, , 20] <- theta_N_init_3[, , 2]
theta_N_init_3[, , 21] <- theta_N_init_3[, , 3]
theta_N_init_3[, , 22] <- theta_N_init_3[, , 4]
theta_N_init_3[, , 23] <- theta_N_init_3[, , 5]
theta_N_init_3[, , 24] <- theta_N_init_3[, , 6]
theta_N_init_3[, , 25] <- theta_N_init_3[, , 7]
theta_N_init_3[, , 26] <- theta_N_init_3[, , 8]
theta_N_init_3[, , 27] <- theta_N_init_3[, , 9]
theta_N_init_3[, , 28] <- theta_N_init_3[, , 10]
theta_N_init_3[, , 29] <- theta_N_init_3[, , 11]
theta_N_init_3[, , 30] <- theta_N_init_3[, , 12]
# -
theta_N_init_3[, , 31] <- theta_N_init_3[, , 13]
theta_N_init_3[, , 32] <- theta_N_init_3[, , 14]
theta_N_init_3[, , 33] <- theta_N_init_3[, , 15]
theta_N_init_3[, , 34] <- theta_N_init_3[, , 16]
theta_N_init_3[, , 35] <- theta_N_init_3[, , 17]
theta_N_init_3[, , 36] <- theta_N_init_3[, , 18]
theta_N_init_3[, , 37] <- theta_N_init_3[, , 1]
theta_N_init_3[, , 38] <- theta_N_init_3[, , 2]
theta_N_init_3[, , 39] <- theta_N_init_3[, , 3]
theta_N_init_3[, , 40] <- theta_N_init_3[, , 4]
theta_N_init_3[, , 41] <- theta_N_init_3[, , 5]
theta_N_init_3[, , 42] <- theta_N_init_3[, , 6]
theta_N_init_3[, , 43] <- theta_N_init_3[, , 7]
theta_N_init_3[, , 44] <- theta_N_init_3[, , 8]
theta_N_init_3[, , 45] <- theta_N_init_3[, , 9]
theta_N_init_3[, , 46] <- theta_N_init_3[, , 10]
theta_N_init_3[, , 47] <- theta_N_init_3[, , 11]
theta_N_init_3[, , 48] <- theta_N_init_3[, , 12]
# Individual claim sizes
theta_X_init_3 <- array(NA, dim = c(nrow(Optim_2$theta_X), K, Dim_init_3))
theta_X_init_3[, , 1] <- cbind(c(0.5 ^ 2 * Optim_2$theta_X[c(1, 2), 1], Optim_2$theta_X[
  -c(1, 2, nrow(Optim_2$theta_X)), 1], 2 * Optim_2$theta_X[nrow(Optim_2$theta_X), 1]), NA, 
  c(2 ^ 2 * Optim_2$theta_X[c(1, 2), 2], Optim_2$theta_X[-c(1, 2, nrow(Optim_2$theta_X)), 2],
    0.5 * Optim_2$theta_N[nrow(Optim_2$theta_N), 2]))
theta_X_init_3[, 2, 1] <- apply(theta_X_init_3[, c(1, 3), 1], 1, mean)
theta_X_init_3[, , 2] <- theta_X_init_3[, c(1, 3, 2), 1]
theta_X_init_3[, , 3] <- theta_X_init_3[, c(2, 1, 3), 1]
theta_X_init_3[, , 4] <- theta_X_init_3[, c(2, 3, 1), 1]
theta_X_init_3[, , 5] <- theta_X_init_3[, c(3, 2, 1), 1]
theta_X_init_3[, , 6] <- theta_X_init_3[, c(3, 1, 2), 1]
theta_X_init_3[, , 7] <- cbind(Optim_2$theta_X[, 1] * c(ifelse(Optim_2$theta_X[
  -nrow(Optim_2$theta_X), 1] < 0, 2, 0.5), 2), NA, Optim_2$theta_X[, 2] * c(ifelse(
    Optim_2$theta_X[-nrow(Optim_2$theta_X), 2] < 0, 0.5, 2), 0.5))
theta_X_init_3[, 2, 7] <- apply(theta_X_init_3[, c(1, 3), 7], 1, mean)
theta_X_init_3[, , 8] <- theta_X_init_3[, c(1, 3, 2), 7]
theta_X_init_3[, , 9] <- theta_X_init_3[, c(2, 1, 3), 7]
theta_X_init_3[, , 10] <- theta_X_init_3[, c(2, 3, 1), 7]
theta_X_init_3[, , 11] <- theta_X_init_3[, c(3, 2, 1), 7]
theta_X_init_3[, , 12] <- theta_X_init_3[, c(3, 1, 2), 7]
# -
theta_X_init_3[, , 13] <- cbind(Optim_2$theta_X[, 1], apply(Optim_2$theta_X, 1, mean),
                                Optim_2$theta_X[, 2])
theta_X_init_3[, , 14] <- theta_X_init_3[, c(1, 3, 2), 13]
theta_X_init_3[, , 15] <- theta_X_init_3[, c(2, 1, 3), 13]
theta_X_init_3[, , 16] <- theta_X_init_3[, c(2, 3, 1), 13]
theta_X_init_3[, , 17] <- theta_X_init_3[, c(3, 2, 1), 13]
theta_X_init_3[, , 18] <- theta_X_init_3[, c(3, 1, 2), 13]
theta_X_init_3[, , 19] <- theta_X_init_3[, , 1]
theta_X_init_3[, , 20] <- theta_X_init_3[, , 2]
theta_X_init_3[, , 21] <- theta_X_init_3[, , 3]
theta_X_init_3[, , 22] <- theta_X_init_3[, , 4]
theta_X_init_3[, , 23] <- theta_X_init_3[, , 5]
theta_X_init_3[, , 24] <- theta_X_init_3[, , 6]
theta_X_init_3[, , 25] <- theta_X_init_3[, , 7]
theta_X_init_3[, , 26] <- theta_X_init_3[, , 8]
theta_X_init_3[, , 27] <- theta_X_init_3[, , 9]
theta_X_init_3[, , 28] <- theta_X_init_3[, , 10]
theta_X_init_3[, , 29] <- theta_X_init_3[, , 11]
theta_X_init_3[, , 30] <- theta_X_init_3[, , 12]
# -
theta_X_init_3[, , 31] <- theta_X_init_3[, , 13]
theta_X_init_3[, , 32] <- theta_X_init_3[, , 14]
theta_X_init_3[, , 33] <- theta_X_init_3[, , 15]
theta_X_init_3[, , 34] <- theta_X_init_3[, , 16]
theta_X_init_3[, , 35] <- theta_X_init_3[, , 17]
theta_X_init_3[, , 36] <- theta_X_init_3[, , 18]
theta_X_init_3[, , 37] <- theta_X_init_3[, , 1]
theta_X_init_3[, , 38] <- theta_X_init_3[, , 2]
theta_X_init_3[, , 39] <- theta_X_init_3[, , 3]
theta_X_init_3[, , 40] <- theta_X_init_3[, , 4]
theta_X_init_3[, , 41] <- theta_X_init_3[, , 5]
theta_X_init_3[, , 42] <- theta_X_init_3[, , 6]
theta_X_init_3[, , 43] <- theta_X_init_3[, , 7]
theta_X_init_3[, , 44] <- theta_X_init_3[, , 8]
theta_X_init_3[, , 45] <- theta_X_init_3[, , 9]
theta_X_init_3[, , 46] <- theta_X_init_3[, , 10]
theta_X_init_3[, , 47] <- theta_X_init_3[, , 11]
theta_X_init_3[, , 48] <- theta_X_init_3[, , 12]
# Ensure constraints are met at initial values
theta_N_init_3[dim(theta_N_init_3)[1], , ] <- 
  pmax(theta_N_init_3[dim(theta_N_init_3)[1], , ], 0.10)
theta_X_init_3[1, , ] <- pmax(theta_X_init_3[1, , ], 0.10)
theta_X_init_3[dim(theta_X_init_3)[1], , ] <- 
  pmax(theta_X_init_3[dim(theta_X_init_3)[1], , ], 1.10)

### Try different starting values
ell_init_3 <- rep(NA, Dim_init_3)
for (init in 1:Dim_init_3) {
  ell_init_3[init] <- tail(EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_3[, , init], 
                                           theta_N_init = theta_N_init_3[, , init],
                                           theta_X_init = theta_X_init_3[, , init], K = K, 
                                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, 
                                           Sizes = MTPL_X_Sev, Ind_i_t = MTPL_Ind_i_t, 
                                           Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW, 
                                           Tol_rel_BW = Tol_rel_BW, Iter_BW = Iter_BW, 
                                           Tol_abs_M = Tol_abs_M, Tol_rel_M = Tol_rel_M, 
                                           Iter_M = Iter_M, Scaled = 1, LB_eps = LB_eps, 
                                           penalty = penalty)$ell_BW, 1)
  print(sprintf(paste('Initial run: %2i / %2i, Log-likelihood: %9.4f'), init, Dim_init_3, 
                ell_init_3[init]))
}
best_init_3 <- which.max(ell_init_3)

### Optimize for K=3 using the Empirical Bayes Baum-Welch algorithm with the BFGS method
Optim_3 <- EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_3[, , best_init_3], 
                           theta_N_init = theta_N_init_3[, , best_init_3],
                           theta_X_init = theta_X_init_3[, , best_init_3], K = K, 
                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, Sizes = MTPL_X_Sev, 
                           Ind_i_t = MTPL_Ind_i_t, Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW_incr, 
                           Tol_rel_BW = Tol_rel_BW_incr, Iter_BW = Iter_BW_incr, 
                           Tol_abs_M = Tol_abs_M_incr, Tol_rel_M = Tol_rel_M_incr, Iter_M = Iter_M, 
                           Scaled = 1, LB_eps = LB_eps, penalty = penalty)

##### Optimize (K=4) #####
### Search grid for initial parameter values for four risk profiles (K=4)
# Number of latent risk profiles
K <- 4
Dim_init_4 <- 36
theta_Z_3_adj <- Optim_3$theta_Z
theta_Z_3_adj[2, ] <- c(0.85, 0.10, 0.05)
## Latent profile assignments
theta_Z_init_4 <- array(NA, dim = c(1 + K, K, Dim_init_4))
theta_Z_init_4[1:K, 1:(K - 1), 1] <- theta_Z_3_adj
theta_Z_init_4[1, , 1] <- c(theta_Z_init_4[1, 1:(K - 1), 1] - c(0.025, 0, 0.025), 0.05)
theta_Z_init_4[2, , 1] <- c(theta_Z_init_4[2, 1:(K - 1), 1] - c(0.05, 0, 0), 0.05)
theta_Z_init_4[3, , 1] <- c(theta_Z_init_4[3, 1:(K - 1), 1] - c(0, 0, 0.05), 0.05)
theta_Z_init_4[4, , 1] <- c(theta_Z_init_4[4, 1:(K - 1), 1] - c(0, 0.025, 0.025), 0.05)
theta_Z_init_4[5, , 1] <- c(0.05, 0.20, 0.30, 0.45)
theta_Z_init_4[, , 2] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 3] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 4] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 5] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 6] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 7] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 8] <- theta_Z_init_4[, , 1]
theta_Z_init_4[, , 9] <- theta_Z_init_4[, , 1]
# -
theta_Z_init_4[1:K, 1:(K - 1), 10] <- theta_Z_3_adj
theta_Z_init_4[1, , 10] <- c(theta_Z_init_4[1, 1:(K - 1), 10] - c(0.05, 0, 0.05), 0.10)
theta_Z_init_4[2, , 10] <- c(theta_Z_init_4[2, 1:(K - 1), 10] - c(0.10, 0, 0), 0.10)
theta_Z_init_4[3, , 10] <- c(theta_Z_init_4[3, 1:(K - 1), 10] - c(0, 0, 0.10), 0.10)
theta_Z_init_4[4, , 10] <- c(theta_Z_init_4[4, 1:(K - 1), 10] - c(0, 0.05, 0.05), 0.10)
theta_Z_init_4[5, , 10] <- c(0.05, 0.20, 0.30, 0.45)
theta_Z_init_4[, , 11] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 12] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 13] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 14] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 15] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 16] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 17] <- theta_Z_init_4[, , 10]
theta_Z_init_4[, , 18] <- theta_Z_init_4[, , 10]
# -
theta_Z_init_4[1:K, 1:(K - 1), 19] <- theta_Z_3_adj
theta_Z_init_4[1, , 19] <- c(theta_Z_init_4[1, 1:(K - 1), 19] - c(0.075, 0.05, 0.075), 0.20)
theta_Z_init_4[2, , 19] <- c(theta_Z_init_4[2, 1:(K - 1), 19] - c(0.20, 0, 0), 0.20)
theta_Z_init_4[3, , 19] <- c(theta_Z_init_4[3, 1:(K - 1), 19] - c(0.025, 0.05, 0.125), 0.20)
theta_Z_init_4[4, , 19] <- c(theta_Z_init_4[4, 1:(K - 1), 19] - c(0, 0.075, 0.125), 0.20)
theta_Z_init_4[5, , 19] <- c(0.05, 0.20, 0.30, 0.45)
theta_Z_init_4[, , 20] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 21] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 22] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 23] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 24] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 25] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 26] <- theta_Z_init_4[, , 19]
theta_Z_init_4[, , 27] <- theta_Z_init_4[, , 19]
# -
theta_Z_init_4[1:K, 1:(K - 1), 28] <- theta_Z_3_adj
theta_Z_init_4[1, , 28] <- c(theta_Z_init_4[1, 1:(K - 1), 28] - c(0.1125, 0.075, 0.1125), 0.30)
theta_Z_init_4[2, , 28] <- c(theta_Z_init_4[2, 1:(K - 1), 28] - c(0.30, 0, 0), 0.30)
theta_Z_init_4[3, , 28] <- c(theta_Z_init_4[3, 1:(K - 1), 28] - c(0.05, 0.075, 0.175), 0.30)
theta_Z_init_4[4, , 28] <- c(theta_Z_init_4[4, 1:(K - 1), 28] - c(0.025, 0.10, 0.175), 0.30)
theta_Z_init_4[5, , 28] <- c(0.05, 0.20, 0.30, 0.45)
theta_Z_init_4[, , 29] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 30] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 31] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 32] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 33] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 34] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 35] <- theta_Z_init_4[, , 28]
theta_Z_init_4[, , 36] <- theta_Z_init_4[, , 28]
## Number of claims
theta_N_init_4 <- array(NA, dim = c(length(Optim_1$theta_N), K, Dim_init_4))
theta_N_init_4[, 1:(K - 1), 1] <- Optim_3$theta_N
theta_N_init_4[, K, 1] <- (Optim_3$theta_N[, 1] + Optim_3$theta_N[, 2]) / 2
theta_N_init_4[, 1:(K - 1), 2] <- Optim_3$theta_N
theta_N_init_4[, K, 2] <- (Optim_3$theta_N[, 2] + Optim_3$theta_N[, 3]) / 2
theta_N_init_4[, 1:(K - 1), 3] <- Optim_3$theta_N
theta_N_init_4[, K, 3] <- 0.85 * Optim_3$theta_N[, 3]
theta_N_init_4[, , 4] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 5] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 6] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 7] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 8] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 9] <- theta_N_init_4[, , 3]
# -
theta_N_init_4[, , 10] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 11] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 12] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 13] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 14] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 15] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 16] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 17] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 18] <- theta_N_init_4[, , 3]
# -
theta_N_init_4[, , 19] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 20] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 21] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 22] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 23] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 24] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 25] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 26] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 27] <- theta_N_init_4[, , 3]
# -
theta_N_init_4[, , 28] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 29] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 30] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 31] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 32] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 33] <- theta_N_init_4[, , 3]
theta_N_init_4[, , 34] <- theta_N_init_4[, , 1]
theta_N_init_4[, , 35] <- theta_N_init_4[, , 2]
theta_N_init_4[, , 36] <- theta_N_init_4[, , 3]
## Individual claim sizes
theta_X_init_4 <- array(NA, dim = c(length(Optim_1$theta_X), K, Dim_init_4))
theta_X_init_4[, 1:(K - 1), 1] <- Optim_3$theta_X
theta_X_init_4[, K, 1] <- (Optim_3$theta_X[, 1] + Optim_3$theta_X[, 2]) / 2
theta_X_init_4[, , 2] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 3] <- theta_X_init_4[, , 1]
theta_X_init_4[, 1:(K - 1), 4] <- Optim_3$theta_X
theta_X_init_4[, K, 4] <- (Optim_3$theta_X[, 2] + Optim_3$theta_X[, 3]) / 2
theta_X_init_4[, , 5] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 6] <- theta_X_init_4[, , 4]
theta_X_init_4[, 1:(K - 1), 7] <- Optim_3$theta_X
theta_X_init_4[, K, 7] <- 1.025 * Optim_3$theta_X[, 3]
theta_X_init_4[, , 8] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 9] <- theta_X_init_4[, , 7]
# -
theta_X_init_4[, , 10] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 11] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 12] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 13] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 14] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 15] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 16] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 17] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 18] <- theta_X_init_4[, , 7]
# -
theta_X_init_4[, , 19] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 20] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 21] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 22] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 23] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 24] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 25] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 26] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 27] <- theta_X_init_4[, , 7]
# -
theta_X_init_4[, , 28] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 29] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 30] <- theta_X_init_4[, , 1]
theta_X_init_4[, , 31] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 32] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 33] <- theta_X_init_4[, , 4]
theta_X_init_4[, , 34] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 35] <- theta_X_init_4[, , 7]
theta_X_init_4[, , 36] <- theta_X_init_4[, , 7]
# Ensure constraints are met at initial values
theta_N_init_4[dim(theta_N_init_4)[1], , ] <- 
  pmax(theta_N_init_4[dim(theta_N_init_4)[1], , ], 0.10)
theta_X_init_4[1, , ] <- pmax(theta_X_init_4[1, , ], 0.10)
theta_X_init_4[dim(theta_X_init_4)[1], , ] <- 
  pmax(theta_X_init_4[dim(theta_X_init_4)[1], , ], 1.10)

### Try different starting values
ell_init_4 <- rep(NA, Dim_init_4)
for (init in 1:Dim_init_4) {
  ell_init_4[init] <- tail(EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_4[, , init], 
                                           theta_N_init = theta_N_init_4[, , init],
                                           theta_X_init = theta_X_init_4[, , init], K = K, 
                                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, 
                                           Sizes = MTPL_X_Sev, Ind_i_t = MTPL_Ind_i_t, 
                                           Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW, 
                                           Tol_rel_BW = Tol_rel_BW, Iter_BW = Iter_BW, 
                                           Tol_abs_M = Tol_abs_M, Tol_rel_M = Tol_rel_M, 
                                           Iter_M = Iter_M, Scaled = 1, LB_eps = LB_eps, 
                                           penalty = penalty)$ell_BW, 1)
  print(sprintf(paste('Initial run: %2i / %2i, Log-likelihood: %9.4f'), init, Dim_init_4, 
                ell_init_4[init]))
}
best_init_4 <- which.max(ell_init_4)

### Additional search grid
Dim_init_4_ext <- 13
theta_Z_init_4_ext <- array(NA, dim = c(1 + K, K, Dim_init_4_ext))
theta_Z_init_4_ext[, , 1] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 1] <- c(0.05, 0.20, 0.30, 0.45)
theta_Z_init_4_ext[, , 2] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 2] <- c(0.05, 0.10, 0.20, 0.65)
theta_Z_init_4_ext[, , 3] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 3] <- c(0.05, 0.20, 0.10, 0.65)
theta_Z_init_4_ext[, , 4] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 4] <- c(0.20, 0.10, 0.10, 0.60)
theta_Z_init_4_ext[, , 5] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 5] <- c(0.40, 0.25, 0.30, 0.15)
theta_Z_init_4_ext[, , 6] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 6] <- c(0.25, 0.40, 0.30, 0.15)
theta_Z_init_4_ext[, , 7] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 7] <- c(0.25, 0.30, 0.40, 0.15)
theta_Z_init_4_ext[, , 8] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 8] <- c(0.40, 0.40, 0.10, 0.10)
theta_Z_init_4_ext[, , 9] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 9] <- c(0.40, 0.10, 0.40, 0.10)
theta_Z_init_4_ext[, , 10] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 10] <- c(0.40, 0.10, 0.10, 0.40)
theta_Z_init_4_ext[, , 11] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 11] <- c(0.10, 0.40, 0.40, 0.10)
theta_Z_init_4_ext[, , 12] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 12] <- c(0.10, 0.40, 0.10, 0.40)
theta_Z_init_4_ext[, , 13] <- theta_Z_init_4[, , best_init_4]
theta_Z_init_4_ext[1 + K, , 13] <- c(0.10, 0.10, 0.40, 0.40)

### Try additional, different starting values
ell_init_4_ext <- rep(NA, Dim_init_4_ext)
for (init in 1:Dim_init_4_ext) {
  ell_init_4_ext[init] <- tail(EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_4_ext[, , init], 
                                               theta_N_init = theta_N_init_4[, , best_init_4],
                                               theta_X_init = theta_X_init_4[, , best_init_4], K = K, 
                                               Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, 
                                               Sizes = MTPL_X_Sev, Ind_i_t = MTPL_Ind_i_t, 
                                               Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                                               RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW, 
                                               Tol_rel_BW = Tol_rel_BW, Iter_BW = Iter_BW, 
                                               Tol_abs_M = Tol_abs_M, Tol_rel_M = Tol_rel_M, 
                                               Iter_M = Iter_M, Scaled = 1, LB_eps = LB_eps, 
                                               penalty = penalty)$ell_BW, 1)
  print(sprintf(paste('Initial run: %2i / %2i, Log-likelihood: %9.4f'), init, Dim_init_4_ext, 
                ell_init_4_ext[init]))
}
best_init_4_ext <- which.max(ell_init_4_ext)

### Optimize for K=4 using the Empirical Bayes Baum-Welch algorithm with the BFGS method
Optim_4 <- EB_BW_BFGS_Full(theta_Z_init = theta_Z_init_4_ext[, , best_init_4_ext], 
                           theta_N_init = theta_N_init_4[, , best_init_4],
                           theta_X_init = theta_X_init_4[, , best_init_4], K = K, 
                           Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, Sizes = MTPL_X_Sev, 
                           Ind_i_t = MTPL_Ind_i_t, Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                           RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW_incr, 
                           Tol_rel_BW = Tol_rel_BW_incr, Iter_BW = Iter_BW_incr, 
                           Tol_abs_M = Tol_abs_M_incr, Tol_rel_M = Tol_rel_M_incr, Iter_M = Iter_M, 
                           Scaled = 1, LB_eps = LB_eps, penalty = penalty)

### Store the final results
