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
### Load results from MTPL_Preprocessing.R

##### Input settings #####
### Load function for Emprical Bayes Baum-Welch algorithm with BFGS method
EB_BW_BFGS_Sparse <- dget('.../Sparse/Functions/EB_BaumWelch_BFGS_Sparse.R')

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

### Remove constant from risk factors for identification purposes
MTPL_RF_Freq <- MTPL_RF_Freq[, -which(colnames(MTPL_RF_Freq) == '(Intercept)')]
MTPL_RF_Sev <- MTPL_RF_Sev[, -which(colnames(MTPL_RF_Sev) == '(Intercept)')]
MTPL_RF_Sev_Freq <- MTPL_RF_Sev_Freq[, -which(colnames(MTPL_RF_Sev_Freq) == '(Intercept)')]

##### Optimize (K=1) #####
### Initial parameter estimates for a single risk profile (K=1)
# Number of latent risk profiles
K <- 1
# Latent profile assignments
theta_Z_init_1 <- NA
# Number of claims
NB_init <- glm.nb(formula = MTPL_N_Freq ~ 1 + MTPL_RF_Freq + offset(log(MTPL_Exp_Freq)), 
                  link = log, control = glm.control(epsilon = Tol_rel_M, maxit = Iter_M))
theta_N_init_1 <- as.vector(c(NB_init$coefficients[-1], NB_init$theta, 
                              NB_init$theta * exp(-NB_init$coefficients[1])))
# Individual claim sizes
G_init <- glm(formula = MTPL_X_Sev ~ 1 + MTPL_RF_Sev, family = Gamma(link = 'log'))
G_init$shape <- gamma.shape(G_init)
GB2_init <- gamlss(formula = MTPL_X_Sev ~ 1, sigma.formula = MTPL_X_Sev ~ 1,
                   nu.formula = MTPL_X_Sev ~ 1 + MTPL_RF_Sev, tau.formula = MTPL_X_Sev ~ 1,
                   family = GB2(mu.link = 'log', sigma.link = 'log',
                                nu.link = 'log', tau.link = 'log'),
                   mu.start = (2 - 1) / G_init$shape$alpha, sigma.start = 1,
                   nu.start = G_init$shape$alpha * as.numeric(G_init$fitted.values),
                   tau.start = 2, sigma.fix = TRUE,
                   control = gamlss.control(c.crit = Tol_abs_M, n.cyc = Iter_M, trace = FALSE), 
                   i.control = glim.control(cc = Tol_abs_M, cyc = Iter_M, bf.cyc = Iter_M,
                                            bf.tol = Tol_abs_M))
theta_X_init_1 <- as.vector(c(exp(GB2_init$nu.coefficients[1]), 
                              GB2_init$nu.coefficients[-1], exp(GB2_init$tau.coefficients),
                              exp(GB2_init$nu.coefficients[1] + GB2_init$mu.coefficients)))

### Optimize for K=1 using the Empirical Bayes Baum-Welch algorithm with the BFGS method
Optim_1 <- EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_1, theta_N_init = theta_N_init_1,
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
Dim_init_2 <- 2 * (8 - 1)
# Latent profile assignments
theta_Z_init_2 <- array(NA, dim = c(1 + K, K, Dim_init_2))
theta_Z_init_2[, , 1] <- t(sapply(1:(K + 1), function(j) rep(1 / K, K)))
theta_Z_init_2[, , 2] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 3] <- theta_Z_init_2[, , 1]
theta_Z_init_2[, , 4] <- theta_Z_init_2[, , 1]
# -
theta_Z_init_2[, , 5] <- rbind(c(0.75, 0.25), c(0.75, 0.25), c(0.25, 0.75))
theta_Z_init_2[, , 6] <- theta_Z_init_2[, , 5]
theta_Z_init_2[, , 7] <- theta_Z_init_2[, , 5]
theta_Z_init_2[, , 8] <- theta_Z_init_2[, , 5]
theta_Z_init_2[, , 9] <- theta_Z_init_2[, , 5]
# -
theta_Z_init_2[, , 10] <- rbind(c(0.85, 0.15), c(0.95, 0.05), c(0.05, 0.95))
theta_Z_init_2[, , 11] <- theta_Z_init_2[, , 10]
theta_Z_init_2[, , 12] <- theta_Z_init_2[, , 10]
theta_Z_init_2[, , 13] <- theta_Z_init_2[, , 10]
theta_Z_init_2[, , 14] <- theta_Z_init_2[, , 10]
# Number of claims
theta_N_init_2 <- matrix(NA, nrow = length(Optim_1$theta_N) + 2 * (K - 1), ncol = Dim_init_2)
theta_N_init_2[, 1] <- c(head(Optim_1$theta_N, -2), 
                         Optim_1$theta_N[length(Optim_1$theta_N) - 1] + K_distort_2 * 
                           sqrt(pmax(Optim_1$V_theta_N[length(Optim_1$theta_N) - 1], Tol_abs_M ^ 2)),
                         Optim_1$theta_N[length(Optim_1$theta_N)] + rev(K_distort_2) * 
                           sqrt(pmax(Optim_1$V_theta_N[length(Optim_1$theta_N)], Tol_abs_M ^ 2)))
theta_N_init_2[, 2] <- theta_N_init_2[, 1]
theta_N_init_2[, 3] <- c(head(Optim_1$theta_N, -2), 
                         c(0.5, 2) * rep(Optim_1$theta_N[length(Optim_1$theta_N) - 1], K),
                         c(2, 0.5) * rep(Optim_1$theta_N[length(Optim_1$theta_N)], K))
theta_N_init_2[, 4] <- theta_N_init_2[, 3]
# -
theta_N_init_2[, 5] <- c(head(Optim_1$theta_N, -2), 
                         rep(Optim_1$theta_N[length(Optim_1$theta_N) - 1], K),
                         rep(Optim_1$theta_N[length(Optim_1$theta_N)], K))
theta_N_init_2[, 6] <- theta_N_init_2[, 1]
theta_N_init_2[, 7] <- theta_N_init_2[, 2]
theta_N_init_2[, 8] <- theta_N_init_2[, 3]
theta_N_init_2[, 9] <- theta_N_init_2[, 4]
# -
theta_N_init_2[, 10] <- theta_N_init_2[, 5]
theta_N_init_2[, 11] <- theta_N_init_2[, 1]
theta_N_init_2[, 12] <- theta_N_init_2[, 2]
theta_N_init_2[, 13] <- theta_N_init_2[, 3]
theta_N_init_2[, 14] <- theta_N_init_2[, 4]
# Individual claim sizes
theta_X_init_2 <- matrix(NA, nrow = length(Optim_1$theta_X) + 2 * (K - 1), ncol = Dim_init_2)
theta_X_init_2[, 1] <- c(head(Optim_1$theta_X, -2), 
                         Optim_1$theta_X[length(Optim_1$theta_X) - 1] + K_distort_2 * 
                           sqrt(pmax(Optim_1$V_theta_X[length(Optim_1$theta_X) - 1], Tol_abs_M ^ 2)),
                         Optim_1$theta_X[length(Optim_1$theta_X)] + rev(K_distort_2) * 
                           sqrt(pmax(Optim_1$V_theta_X[length(Optim_1$theta_X)], Tol_abs_M ^ 2)))
theta_X_init_2[, 2] <- theta_X_init_2[c(1:(nrow(theta_X_init_2) - 2 * K), nrow(theta_X_init_2) - K, 
                                        nrow(theta_X_init_2) - K - 1, nrow(theta_X_init_2), 
                                        nrow(theta_X_init_2) - K + 1), 1]
theta_X_init_2[, 3] <- c(head(Optim_1$theta_X, -2), 
                         c(2, 0.5) * rep(Optim_1$theta_X[length(Optim_1$theta_X) - 1], K),
                         c(0.5, 2) * rep(Optim_1$theta_X[length(Optim_1$theta_X)], K))
theta_X_init_2[, 4] <- theta_X_init_2[c(1:(nrow(theta_X_init_2) - 2 * K), nrow(theta_X_init_2) - K, 
                                        nrow(theta_X_init_2) - K - 1, nrow(theta_X_init_2), 
                                        nrow(theta_X_init_2) - K + 1), 3]
# -
theta_X_init_2[, 5] <- c(head(Optim_1$theta_X, -2), 
                         rep(Optim_1$theta_X[length(Optim_1$theta_X) - 1], K),
                         rep(Optim_1$theta_X[length(Optim_1$theta_X)], K))
theta_X_init_2[, 6] <- theta_X_init_2[, 1]
theta_X_init_2[, 7] <- theta_X_init_2[, 2]
theta_X_init_2[, 8] <- theta_X_init_2[, 3]
theta_X_init_2[, 9] <- theta_X_init_2[, 4]
# -
theta_X_init_2[, 10] <- theta_X_init_2[, 5]
theta_X_init_2[, 11] <- theta_X_init_2[, 1]
theta_X_init_2[, 12] <- theta_X_init_2[, 2]
theta_X_init_2[, 13] <- theta_X_init_2[, 3]
theta_X_init_2[, 14] <- theta_X_init_2[, 4]
# Ensure constraints are met at initial values
theta_N_init_2[(nrow(theta_N_init_2) - 2 * K + 1):nrow(theta_N_init_2), ] <- 
  pmax(theta_N_init_2[(nrow(theta_N_init_2) - 2 * K + 1):nrow(theta_N_init_2), ], 0.10)
theta_X_init_2[(nrow(theta_X_init_2) - 2 * K + 1):(nrow(theta_X_init_2) - K), ] <- 
  pmax(theta_X_init_2[(nrow(theta_X_init_2) - 2 * K + 1):(nrow(theta_X_init_2) - K), ], 1.10)
theta_X_init_2[(nrow(theta_X_init_2) - K + 1):nrow(theta_X_init_2), ] <- 
  pmax(theta_X_init_2[(nrow(theta_X_init_2) - K + 1):nrow(theta_X_init_2), ], 0.10)

### Try different starting values
ell_init_2 <- rep(NA, Dim_init_2)
for (init in 1:Dim_init_2) {
  ell_init_2[init] <- tail(EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_2[, , init], 
                                             theta_N_init = theta_N_init_2[, init],
                                             theta_X_init = theta_X_init_2[, init], K = K, 
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
Optim_2 <- EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_2[, , best_init_2], 
                             theta_N_init = theta_N_init_2[, best_init_2],
                             theta_X_init = theta_X_init_2[, best_init_2], K = K, 
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
# Latent profile assignments
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
# Number of claims
theta_N_init_3 <- matrix(NA, nrow = length(Optim_1$theta_N) + 2 * (K - 1), ncol = Dim_init_3)
theta_N_init_3[, 1] <- c(head(Optim_2$theta_N, -4), tail(Optim_2$theta_N, 4)[c(1, 1, 2)] + 
                           c(K_distort_3[1], NA, K_distort_3[3]) * 
                           sqrt(pmax(tail(Optim_2$V_theta_N, 4)[c(1, 1, 2)], Tol_abs_M ^ 2)),
                         tail(Optim_2$theta_N, 4)[c(3, 3, 4)] + 
                           c(K_distort_3[3], NA, K_distort_3[1]) * 
                           sqrt(pmax(tail(Optim_2$V_theta_N, 4)[c(3, 3, 4)], Tol_abs_M ^ 2)))
theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 2, nrow(theta_N_init_3) - 1), 1] <- 
  c(mean(theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 1, nrow(theta_N_init_3) - K), 1]), 
    mean(theta_N_init_3[c(nrow(theta_N_init_3) - K + 1, nrow(theta_N_init_3)), 1]))
theta_N_init_3[, 2] <- theta_N_init_3[, 1]
theta_N_init_3[, 3] <- theta_N_init_3[, 2]
theta_N_init_3[, 4] <- theta_N_init_3[, 3]
theta_N_init_3[, 5] <- theta_N_init_3[, 4]
theta_N_init_3[, 6] <- theta_N_init_3[, 5]
theta_N_init_3[, 7] <- c(head(Optim_2$theta_N, -4), c(0.5, NA, 2, 2, NA, 0.5) * 
                           tail(Optim_2$theta_N, 4)[c(1, 1, 2, 3, 3, 4)])
theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 2, nrow(theta_N_init_3) - 1), 7] <- 
  c(mean(theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 1, nrow(theta_N_init_3) - K), 7]),
    mean(theta_N_init_3[c(nrow(theta_N_init_3) - K + 1, nrow(theta_N_init_3)), 7]))
theta_N_init_3[, 8] <- theta_N_init_3[, 7]
theta_N_init_3[, 9] <- theta_N_init_3[, 8]
theta_N_init_3[, 10] <- theta_N_init_3[, 9]
theta_N_init_3[, 11] <- theta_N_init_3[, 10]
theta_N_init_3[, 12] <- theta_N_init_3[, 11]
# -
theta_N_init_3[, 13] <- c(head(Optim_2$theta_N, -4), tail(Optim_2$theta_N, 4)[c(1, 1, 2, 3, 3, 4)])
theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 2, nrow(theta_N_init_3) - 1), 13] <- 
  c(mean(theta_N_init_3[c(nrow(theta_N_init_3) - 2 * K + 1, nrow(theta_N_init_3) - K), 13]), 
    mean(theta_N_init_3[c(nrow(theta_N_init_3) - K + 1, nrow(theta_N_init_3)), 13]))
theta_N_init_3[, 14] <- theta_N_init_3[, 13]
theta_N_init_3[, 15] <- theta_N_init_3[, 14]
theta_N_init_3[, 16] <- theta_N_init_3[, 15]
theta_N_init_3[, 17] <- theta_N_init_3[, 16]
theta_N_init_3[, 18] <- theta_N_init_3[, 17]
theta_N_init_3[, 19] <- theta_N_init_3[, 1]
theta_N_init_3[, 20] <- theta_N_init_3[, 2]
theta_N_init_3[, 21] <- theta_N_init_3[, 3]
theta_N_init_3[, 22] <- theta_N_init_3[, 4]
theta_N_init_3[, 23] <- theta_N_init_3[, 5]
theta_N_init_3[, 24] <- theta_N_init_3[, 6]
theta_N_init_3[, 25] <- theta_N_init_3[, 7]
theta_N_init_3[, 26] <- theta_N_init_3[, 8]
theta_N_init_3[, 27] <- theta_N_init_3[, 9]
theta_N_init_3[, 28] <- theta_N_init_3[, 10]
theta_N_init_3[, 29] <- theta_N_init_3[, 11]
theta_N_init_3[, 30] <- theta_N_init_3[, 12]
# -
theta_N_init_3[, 31] <- theta_N_init_3[, 13]
theta_N_init_3[, 32] <- theta_N_init_3[, 14]
theta_N_init_3[, 33] <- theta_N_init_3[, 15]
theta_N_init_3[, 34] <- theta_N_init_3[, 16]
theta_N_init_3[, 35] <- theta_N_init_3[, 17]
theta_N_init_3[, 36] <- theta_N_init_3[, 18]
theta_N_init_3[, 37] <- theta_N_init_3[, 1]
theta_N_init_3[, 38] <- theta_N_init_3[, 2]
theta_N_init_3[, 39] <- theta_N_init_3[, 3]
theta_N_init_3[, 40] <- theta_N_init_3[, 4]
theta_N_init_3[, 41] <- theta_N_init_3[, 5]
theta_N_init_3[, 42] <- theta_N_init_3[, 6]
theta_N_init_3[, 43] <- theta_N_init_3[, 7]
theta_N_init_3[, 44] <- theta_N_init_3[, 8]
theta_N_init_3[, 45] <- theta_N_init_3[, 9]
theta_N_init_3[, 46] <- theta_N_init_3[, 10]
theta_N_init_3[, 47] <- theta_N_init_3[, 11]
theta_N_init_3[, 48] <- theta_N_init_3[, 12]
# Individual claim sizes
theta_X_init_3 <- matrix(NA, nrow = length(Optim_1$theta_X) + 2 * (K - 1), ncol = Dim_init_3)
theta_X_init_3[, 1] <- c(head(Optim_2$theta_X, -4), tail(Optim_2$theta_X, 4)[c(1, 1, 2)] + 
                           c(K_distort_3[3], NA, K_distort_3[1]) * 
                           sqrt(pmax(tail(Optim_2$V_theta_X, 4)[c(1, 1, 2)], Tol_abs_M ^ 2)),
                         tail(Optim_2$theta_X, 4)[c(3, 3, 4)] + 
                           c(K_distort_3[1], NA, K_distort_3[3]) * 
                           sqrt(pmax(tail(Optim_2$V_theta_X, 4)[c(3, 3, 4)], Tol_abs_M ^ 2)))
theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 2, nrow(theta_X_init_3) - 1), 1] <- 
  c(mean(theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 1, nrow(theta_X_init_3) - K), 1]), 
    mean(theta_X_init_3[c(nrow(theta_X_init_3) - K + 1, nrow(theta_X_init_3)), 1]))
theta_X_init_3[, 2] <- c(head(theta_X_init_3[, 1], - 2 * K), 
                         tail(theta_X_init_3[, 1], 2 * K)[c(1, 3, 2, 3 + c(1, 3, 2))])
theta_X_init_3[, 3] <- c(head(theta_X_init_3[, 1], - 2 * K), 
                         tail(theta_X_init_3[, 1], 2 * K)[c(2, 1, 3, 3 + c(2, 1, 3))])
theta_X_init_3[, 4] <- c(head(theta_X_init_3[, 1], - 2 * K), 
                         tail(theta_X_init_3[, 1], 2 * K)[c(2, 3, 1, 3 + c(2, 3, 1))])
theta_X_init_3[, 5] <- c(head(theta_X_init_3[, 1], - 2 * K), 
                         tail(theta_X_init_3[, 1], 2 * K)[c(3, 2, 1, 3 + c(3, 2, 1))])
theta_X_init_3[, 6] <- c(head(theta_X_init_3[, 1], - 2 * K), 
                         tail(theta_X_init_3[, 1], 2 * K)[c(3, 1, 2, 3 + c(3, 1, 2))])
theta_X_init_3[, 7] <- c(head(Optim_2$theta_X, -4), c(2, NA, 0.5, 0.5, NA, 2) * 
                           tail(Optim_2$theta_X, 4)[c(1, 1, 2, 3, 3, 4)])
theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 2, nrow(theta_X_init_3) - 1), 7] <- 
  c(mean(theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 1, nrow(theta_X_init_3) - K), 7]), 
    mean(theta_X_init_3[c(nrow(theta_X_init_3) - K + 1, nrow(theta_X_init_3)), 7]))
theta_X_init_3[, 8] <- c(head(theta_X_init_3[, 7], - 2 * K), 
                         tail(theta_X_init_3[, 7], 2 * K)[c(1, 3, 2, 3 + c(1, 3, 2))])
theta_X_init_3[, 9] <- c(head(theta_X_init_3[, 7], - 2 * K), 
                         tail(theta_X_init_3[, 7], 2 * K)[c(2, 1, 3, 3 + c(2, 1, 3))])
theta_X_init_3[, 10] <- c(head(theta_X_init_3[, 7], - 2 * K), 
                          tail(theta_X_init_3[, 7], 2 * K)[c(2, 3, 1, 3 + c(2, 3, 1))])
theta_X_init_3[, 11] <- c(head(theta_X_init_3[, 7], - 2 * K), 
                          tail(theta_X_init_3[, 7], 2 * K)[c(3, 2, 1, 3 + c(3, 2, 1))])
theta_X_init_3[, 12] <- c(head(theta_X_init_3[, 7], - 2 * K), 
                          tail(theta_X_init_3[, 7], 2 * K)[c(3, 1, 2, 3 + c(3, 1, 2))])
# -
theta_X_init_3[, 13] <- c(head(Optim_2$theta_X, -4), tail(Optim_2$theta_X, 4)[c(1, 1, 2, 3, 3, 4)])
theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 2, nrow(theta_X_init_3) - 1), 13] <- 
  c(mean(theta_X_init_3[c(nrow(theta_X_init_3) - 2 * K + 1, nrow(theta_X_init_3) - K), 13]), 
    mean(theta_X_init_3[c(nrow(theta_X_init_3) - K + 1, nrow(theta_X_init_3)), 13]))
theta_X_init_3[, 14] <- c(head(theta_X_init_3[, 13], - 2 * K), 
                          tail(theta_X_init_3[, 13], 2 * K)[c(1, 3, 2, 3 + c(1, 3, 2))])
theta_X_init_3[, 15] <- c(head(theta_X_init_3[, 13], - 2 * K), 
                          tail(theta_X_init_3[, 13], 2 * K)[c(2, 1, 3, 3 + c(2, 1, 3))])
theta_X_init_3[, 16] <- c(head(theta_X_init_3[, 13], - 2 * K), 
                          tail(theta_X_init_3[, 13], 2 * K)[c(2, 3, 1, 3 + c(2, 3, 1))])
theta_X_init_3[, 17] <- c(head(theta_X_init_3[, 13], - 2 * K), 
                          tail(theta_X_init_3[, 13], 2 * K)[c(3, 2, 1, 3 + c(3, 2, 1))])
theta_X_init_3[, 18] <- c(head(theta_X_init_3[, 13], - 2 * K), 
                          tail(theta_X_init_3[, 13], 2 * K)[c(3, 1, 2, 3 + c(3, 1, 2))])
theta_X_init_3[, 19] <- theta_X_init_3[, 1]
theta_X_init_3[, 20] <- theta_X_init_3[, 2]
theta_X_init_3[, 21] <- theta_X_init_3[, 3]
theta_X_init_3[, 22] <- theta_X_init_3[, 4]
theta_X_init_3[, 23] <- theta_X_init_3[, 5]
theta_X_init_3[, 24] <- theta_X_init_3[, 6]
theta_X_init_3[, 25] <- theta_X_init_3[, 7]
theta_X_init_3[, 26] <- theta_X_init_3[, 8]
theta_X_init_3[, 27] <- theta_X_init_3[, 9]
theta_X_init_3[, 28] <- theta_X_init_3[, 10]
theta_X_init_3[, 29] <- theta_X_init_3[, 11]
theta_X_init_3[, 30] <- theta_X_init_3[, 12]
# -
theta_X_init_3[, 31] <- theta_X_init_3[, 13]
theta_X_init_3[, 32] <- theta_X_init_3[, 14]
theta_X_init_3[, 33] <- theta_X_init_3[, 15]
theta_X_init_3[, 34] <- theta_X_init_3[, 16]
theta_X_init_3[, 35] <- theta_X_init_3[, 17]
theta_X_init_3[, 36] <- theta_X_init_3[, 18]
theta_X_init_3[, 37] <- theta_X_init_3[, 1]
theta_X_init_3[, 38] <- theta_X_init_3[, 2]
theta_X_init_3[, 39] <- theta_X_init_3[, 3]
theta_X_init_3[, 40] <- theta_X_init_3[, 4]
theta_X_init_3[, 41] <- theta_X_init_3[, 5]
theta_X_init_3[, 42] <- theta_X_init_3[, 6]
theta_X_init_3[, 43] <- theta_X_init_3[, 7]
theta_X_init_3[, 44] <- theta_X_init_3[, 8]
theta_X_init_3[, 45] <- theta_X_init_3[, 9]
theta_X_init_3[, 46] <- theta_X_init_3[, 10]
theta_X_init_3[, 47] <- theta_X_init_3[, 11]
theta_X_init_3[, 48] <- theta_X_init_3[, 12]
# Ensure constraints are met at initial values
theta_N_init_3[(nrow(theta_N_init_3) - 2 * K + 1):nrow(theta_N_init_3), ] <- 
  pmax(theta_N_init_3[(nrow(theta_N_init_3) - 2 * K + 1):nrow(theta_N_init_3), ], 0.10)
theta_X_init_3[(nrow(theta_X_init_3) - 2 * K + 1):(nrow(theta_X_init_3) - K), ] <- 
  pmax(theta_X_init_3[(nrow(theta_X_init_3) - 2 * K + 1):(nrow(theta_X_init_3) - K), ], 1.10)
theta_X_init_3[(nrow(theta_X_init_3) - K + 1):nrow(theta_X_init_3), ] <- 
  pmax(theta_X_init_3[(nrow(theta_X_init_3) - K + 1):nrow(theta_X_init_3), ], 0.10)

### Try different starting values
ell_init_3 <- rep(NA, Dim_init_3)
for (init in 1:Dim_init_3) {
  ell_init_3[init] <- tail(EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_3[, , init], 
                                             theta_N_init = theta_N_init_3[, init],
                                             theta_X_init = theta_X_init_3[, init], K = K, 
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
Optim_3 <- EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_3[, , best_init_3], 
                             theta_N_init = theta_N_init_3[, best_init_3],
                             theta_X_init = theta_X_init_3[, best_init_3], K = K, 
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
## Latent profile assignments
theta_Z_init_4 <- array(NA, dim = c(1 + K, K, Dim_init_4))
theta_Z_init_4[1:K, 1:(K - 1), 1] <- theta_Z_3_adj
theta_Z_init_4[1, , 1] <- c(theta_Z_init_4[1, 1:(K - 1), 1] - c(0.05, 0, 0), 0.05)
theta_Z_init_4[2, , 1] <- c(theta_Z_init_4[2, 1:(K - 1), 1] - c(0.05, 0, 0), 0.05)
theta_Z_init_4[3, , 1] <- c(theta_Z_init_4[3, 1:(K - 1), 1] - c(0, 0.05, 0), 0.05)
theta_Z_init_4[4, , 1] <- c(theta_Z_init_4[4, 1:(K - 1), 1] - c(0, 0, 0.05), 0.05)
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
theta_Z_init_4[1, , 10] <- c(theta_Z_init_4[1, 1:(K - 1), 10] - c(0.10, 0, 0), 0.10)
theta_Z_init_4[2, , 10] <- c(theta_Z_init_4[2, 1:(K - 1), 10] - c(0.10, 0, 0), 0.10)
theta_Z_init_4[3, , 10] <- c(theta_Z_init_4[3, 1:(K - 1), 10] - c(0, 0.10, 0), 0.10)
theta_Z_init_4[4, , 10] <- c(theta_Z_init_4[4, 1:(K - 1), 10] - c(0, 0, 0.10), 0.10)
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
theta_Z_init_4[1, , 19] <- c(theta_Z_init_4[1, 1:(K - 1), 19] - c(0.15, 0.025, 0.025), 0.20)
theta_Z_init_4[2, , 19] <- c(theta_Z_init_4[2, 1:(K - 1), 19] - c(0.20, 0, 0), 0.20)
theta_Z_init_4[3, , 19] <- c(theta_Z_init_4[3, 1:(K - 1), 19] - c(0, 0.20, 0), 0.20)
theta_Z_init_4[4, , 19] <- c(theta_Z_init_4[4, 1:(K - 1), 19] - c(0, 0, 0.20), 0.20)
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
theta_Z_init_4[1, , 28] <- c(theta_Z_init_4[1, 1:(K - 1), 28] - c(0.15, 0.075, 0.075), 0.30)
theta_Z_init_4[2, , 28] <- c(theta_Z_init_4[2, 1:(K - 1), 28] - c(0.30, 0, 0), 0.30)
theta_Z_init_4[3, , 28] <- c(theta_Z_init_4[3, 1:(K - 1), 28] - c(0.025, 0.25, 0.025), 0.30)
theta_Z_init_4[4, , 28] <- c(theta_Z_init_4[4, 1:(K - 1), 28] - c(0.025, 0.025, 0.25), 0.30)
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
theta_N_init_4 <- matrix(NA, nrow = length(Optim_1$theta_N) + 2 * (K - 1), ncol = Dim_init_4)
theta_N_init_4[, 1] <- c(head(Optim_3$theta_N, -(K - 1)), NA, tail(Optim_3$theta_N, K - 1), NA)
theta_N_init_4[c(length(Optim_1$theta_N) - 1 + K - 1, nrow(theta_N_init_4)), 1] <- 
  c(Optim_3$theta_N[length(Optim_1$theta_N) - 1] + Optim_3$theta_N[length(Optim_1$theta_N)],
    Optim_3$theta_N[length(Optim_1$theta_N) + 2] + Optim_3$theta_N[length(Optim_1$theta_N) + 3]) / 2
theta_N_init_4[, 2] <- c(head(Optim_3$theta_N, -(K - 1)), NA, tail(Optim_3$theta_N, K - 1), NA)
theta_N_init_4[c(length(Optim_1$theta_N) - 1 + K - 1, nrow(theta_N_init_4)), 2] <- 
  c(Optim_3$theta_N[length(Optim_1$theta_N)] + Optim_3$theta_N[length(Optim_1$theta_N) + 1],
    Optim_3$theta_N[length(Optim_1$theta_N) + 3] + Optim_3$theta_N[length(Optim_1$theta_N) + 4]) / 2
theta_N_init_4[, 3] <- c(head(Optim_3$theta_N, -(K - 1)), NA, tail(Optim_3$theta_N, K - 1), NA)
theta_N_init_4[c(length(Optim_1$theta_N) - 1 + K - 1, nrow(theta_N_init_4)), 3] <- 
  c(sqrt(1.35) * Optim_3$theta_N[length(Optim_1$theta_N) + 1],
    1 / sqrt(1.35) * Optim_3$theta_N[length(Optim_1$theta_N) + 4])
theta_N_init_4[, 4] <- theta_N_init_4[, 1]
theta_N_init_4[, 5] <- theta_N_init_4[, 2]
theta_N_init_4[, 6] <- theta_N_init_4[, 3]
theta_N_init_4[, 7] <- theta_N_init_4[, 1]
theta_N_init_4[, 8] <- theta_N_init_4[, 2]
theta_N_init_4[, 9] <- theta_N_init_4[, 3]
# -
theta_N_init_4[, 10] <- theta_N_init_4[, 1]
theta_N_init_4[, 11] <- theta_N_init_4[, 2]
theta_N_init_4[, 12] <- theta_N_init_4[, 3]
theta_N_init_4[, 13] <- theta_N_init_4[, 1]
theta_N_init_4[, 14] <- theta_N_init_4[, 2]
theta_N_init_4[, 15] <- theta_N_init_4[, 3]
theta_N_init_4[, 16] <- theta_N_init_4[, 1]
theta_N_init_4[, 17] <- theta_N_init_4[, 2]
theta_N_init_4[, 18] <- theta_N_init_4[, 3]
# -
theta_N_init_4[, 19] <- theta_N_init_4[, 1]
theta_N_init_4[, 20] <- theta_N_init_4[, 2]
theta_N_init_4[, 21] <- theta_N_init_4[, 3]
theta_N_init_4[, 22] <- theta_N_init_4[, 1]
theta_N_init_4[, 23] <- theta_N_init_4[, 2]
theta_N_init_4[, 24] <- theta_N_init_4[, 3]
theta_N_init_4[, 25] <- theta_N_init_4[, 1]
theta_N_init_4[, 26] <- theta_N_init_4[, 2]
theta_N_init_4[, 27] <- theta_N_init_4[, 3]
# -
theta_N_init_4[, 28] <- theta_N_init_4[, 1]
theta_N_init_4[, 29] <- theta_N_init_4[, 2]
theta_N_init_4[, 30] <- theta_N_init_4[, 3]
theta_N_init_4[, 31] <- theta_N_init_4[, 1]
theta_N_init_4[, 32] <- theta_N_init_4[, 2]
theta_N_init_4[, 33] <- theta_N_init_4[, 3]
theta_N_init_4[, 34] <- theta_N_init_4[, 1]
theta_N_init_4[, 35] <- theta_N_init_4[, 2]
theta_N_init_4[, 36] <- theta_N_init_4[, 3]
## Individual claim sizes
theta_X_init_4 <- matrix(NA, nrow = length(Optim_1$theta_X) + 2 * (K - 1), ncol = Dim_init_4)
theta_X_init_4[, 1] <- c(head(Optim_3$theta_X, -(K - 1)), NA, tail(Optim_3$theta_X, K - 1), NA)
theta_X_init_4[c(length(Optim_1$theta_X) - 1 + K - 1, nrow(theta_X_init_4)), 1] <- 
  c(Optim_3$theta_X[length(Optim_1$theta_X) - 1] + Optim_3$theta_X[length(Optim_1$theta_X)],
    Optim_3$theta_X[length(Optim_1$theta_X) + 2] + Optim_3$theta_X[length(Optim_1$theta_X) + 3]) / 2
theta_X_init_4[, 2] <- theta_X_init_4[, 1]
theta_X_init_4[, 3] <- theta_X_init_4[, 1]
theta_X_init_4[, 4] <- c(head(Optim_3$theta_X, -(K - 1)), NA, tail(Optim_3$theta_X, K - 1), NA)
theta_X_init_4[c(length(Optim_1$theta_X) - 1 + K - 1, nrow(theta_X_init_4)), 4] <- 
  c(Optim_3$theta_X[length(Optim_1$theta_X)] + Optim_3$theta_X[length(Optim_1$theta_X) + 1],
    Optim_3$theta_X[length(Optim_1$theta_X) + 3] + Optim_3$theta_X[length(Optim_1$theta_X) + 4]) / 2
theta_X_init_4[, 5] <- theta_X_init_4[, 2]
theta_X_init_4[, 6] <- theta_X_init_4[, 2]
theta_X_init_4[, 7] <- c(head(Optim_3$theta_X, -(K - 1)), NA, tail(Optim_3$theta_X, K - 1), NA)
theta_X_init_4[c(length(Optim_1$theta_X) - 1 + K - 1, nrow(theta_X_init_4)), 7] <- 
  c(sqrt(1.25) * Optim_3$theta_X[length(Optim_1$theta_X) + 1], 
    1 / sqrt(1.25) * Optim_3$theta_X[length(Optim_1$theta_X) + 4])
theta_X_init_4[, 8] <- theta_X_init_4[, 3]
theta_X_init_4[, 9] <- theta_X_init_4[, 3]
# -
theta_X_init_4[, 10] <- theta_X_init_4[, 1]
theta_X_init_4[, 11] <- theta_X_init_4[, 1]
theta_X_init_4[, 12] <- theta_X_init_4[, 1]
theta_X_init_4[, 13] <- theta_X_init_4[, 4]
theta_X_init_4[, 14] <- theta_X_init_4[, 4]
theta_X_init_4[, 15] <- theta_X_init_4[, 4]
theta_X_init_4[, 16] <- theta_X_init_4[, 7]
theta_X_init_4[, 17] <- theta_X_init_4[, 7]
theta_X_init_4[, 18] <- theta_X_init_4[, 7]
# -
theta_X_init_4[, 19] <- theta_X_init_4[, 1]
theta_X_init_4[, 20] <- theta_X_init_4[, 1]
theta_X_init_4[, 21] <- theta_X_init_4[, 1]
theta_X_init_4[, 22] <- theta_X_init_4[, 4]
theta_X_init_4[, 23] <- theta_X_init_4[, 4]
theta_X_init_4[, 24] <- theta_X_init_4[, 4]
theta_X_init_4[, 25] <- theta_X_init_4[, 7]
theta_X_init_4[, 26] <- theta_X_init_4[, 7]
theta_X_init_4[, 27] <- theta_X_init_4[, 7]
# -
theta_X_init_4[, 28] <- theta_X_init_4[, 1]
theta_X_init_4[, 29] <- theta_X_init_4[, 1]
theta_X_init_4[, 30] <- theta_X_init_4[, 1]
theta_X_init_4[, 31] <- theta_X_init_4[, 4]
theta_X_init_4[, 32] <- theta_X_init_4[, 4]
theta_X_init_4[, 33] <- theta_X_init_4[, 4]
theta_X_init_4[, 34] <- theta_X_init_4[, 7]
theta_X_init_4[, 35] <- theta_X_init_4[, 7]
theta_X_init_4[, 36] <- theta_X_init_4[, 7]
# Ensure constraints are met at initial values
theta_N_init_4[(nrow(theta_N_init_4) - 2 * K + 1):nrow(theta_N_init_4), ] <- 
  pmax(theta_N_init_4[(nrow(theta_N_init_4) - 2 * K + 1):nrow(theta_N_init_4), ], 0.10)
theta_X_init_4[(nrow(theta_X_init_4) - 2 * K + 1):(nrow(theta_X_init_4) - K), ] <- 
  pmax(theta_X_init_4[(nrow(theta_X_init_4) - 2 * K + 1):(nrow(theta_X_init_4) - K), ], 1.10)
theta_X_init_4[(nrow(theta_X_init_4) - K + 1):nrow(theta_X_init_4), ] <- 
  pmax(theta_X_init_4[(nrow(theta_X_init_4) - K + 1):nrow(theta_X_init_4), ], 0.10)

### Try different starting values
ell_init_4 <- rep(NA, Dim_init_4)
for (init in 1:Dim_init_4) {
  ell_init_4[init] <- tail(EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_4[, , init], 
                                             theta_N_init = theta_N_init_4[, init],
                                             theta_X_init = theta_X_init_4[, init], K = K, 
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
  ell_init_4_ext[init] <- tail(EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_4_ext[, , init], 
                                                 theta_N_init = theta_N_init_4[, best_init_4],
                                                 theta_X_init = theta_X_init_4[, best_init_4], K = K, 
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
Optim_4 <- EB_BW_BFGS_Sparse(theta_Z_init = theta_Z_init_4_ext[, , best_init_4_ext], 
                             theta_N_init = theta_N_init_4[, best_init_4],
                             theta_X_init = theta_X_init_4[, best_init_4], K = K, 
                             Counts = MTPL_N_Freq, Expo = MTPL_Exp_Freq, Sizes = MTPL_X_Sev, 
                             Ind_i_t = MTPL_Ind_i_t, Ind_t_Ti = MTPL_Ind_t_Ti, RF_N = MTPL_RF_Freq, 
                             RF_X = MTPL_RF_Sev, Tol_abs_BW = Tol_abs_BW_incr, 
                             Tol_rel_BW = Tol_rel_BW_incr, Iter_BW = Iter_BW_incr, 
                             Tol_abs_M = Tol_abs_M_incr, Tol_rel_M = Tol_rel_M_incr, Iter_M = Iter_M, 
                             Scaled = 1, LB_eps = LB_eps, penalty = penalty)

### Store the final results
