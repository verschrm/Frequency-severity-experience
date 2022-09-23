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
Packages <- c('dplyr', 'GB2')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load results from MTPL_Bayesian_HMM_Full.R

##### Predictions K=1 #####
### Estimated prior parameters
cbind(tail(Optim_1$theta_N, 1), tail(Optim_1$theta_X, 1))

### Estimated prior transition probabilities
100 * Optim_1$theta_Z

### Mean predictors
# Prior mean predictor
lamb_1 <- pmin(pmax(MTPL_Exp_Freq * exp(MTPL_RF_Freq %*% Optim_1$theta_N[
  1:ncol(MTPL_RF_Freq)]), LB_eps), 1e+15) / MTPL_Exp_Freq
mu_phi_1 <- pmin(pmax(Optim_1$theta_X[1] * exp(MTPL_RF_Sev_Freq %*% Optim_1$theta_X[
  2:(1 + ncol(MTPL_RF_Sev))]), LB_eps), 1e+15) / Optim_1$theta_X[1]
# Posterior parameters
expo_lamb_1 <- pmin(pmax(MTPL_Exp_Freq * lamb_1, LB_eps), 1e+15)
mu_1 <- pmin(pmax(Optim_1$theta_X[1] * mu_phi_1, LB_eps), 1e+15)
Update_alpha_U_1 <- rep(0, nrow(MTPL_RF_Freq))
Update_beta_U_1 <- rep(0, nrow(MTPL_RF_Freq))
Update_alpha_V_1 <- rep(0, nrow(MTPL_RF_Freq))
Update_beta_V_1 <- rep(0, nrow(MTPL_RF_Freq))
alpha_N_1 <- rep(0, nrow(MTPL_RF_Freq))
alpha_N_1[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  mu_1[MTPL_Ind_i_t] ~ MTPL_Ind_i_t, FUN = sum)[, -1])
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
for (t in 2:max(MTPL_Ind_t_Ti[, 2])){
  t_Current <- which(MTPL_Ind_t_Ti[, 1] == t)
  t_Previous <- which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t))
  Update_alpha_U_1[t_Current] <- Update_alpha_U_1[t_Previous] + MTPL_N_Freq[t_Previous]
  Update_beta_U_1[t_Current] <- Update_beta_U_1[t_Previous] + expo_lamb_1[t_Previous]
  Update_alpha_V_1[t_Current] <- Update_alpha_V_1[t_Previous] + alpha_N_1[t_Previous]
  Update_beta_V_1[t_Current] <- Update_beta_V_1[t_Previous] + Sizes_N[t_Previous] * 
    Optim_1$theta_X[1]
}
# Posterior mean predictor
lamb_post_1 <- lamb_1 * (tail(Optim_1$theta_N, 1) + Update_alpha_U_1) /
  (tail(Optim_1$theta_N, 1) + Update_beta_U_1)
mu_phi_post_1 <- mu_phi_1 * (tail(Optim_1$theta_X, 1) - 1 + Update_beta_V_1) /
  (tail(Optim_1$theta_X, 1) - 1 + Update_alpha_V_1)

### Response predictions
# Mean claim frequency
rbind(c(sum(MTPL_N_Freq), MTPL_Exp_Freq %*% lamb_1, MTPL_Exp_Freq %*% lamb_post_2) / 
        sum(MTPL_Exp_Freq),
      c(weightedMedian(x = MTPL_N_Freq / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = lamb_2, w = MTPL_Exp_Freq),
        weightedMedian(x = lamb_post_2, w = MTPL_Exp_Freq)))
# Mean [median] claim severity
rbind(c(mean(MTPL_X_Sev), mean(mu_phi_1[MTPL_Ind_i_t, ]), 
        mean(mu_phi_post_1[MTPL_Ind_i_t, ])), 
      c(median(MTPL_X_Sev), median(mu_phi_1[MTPL_Ind_i_t, ]), 
        median(mu_phi_post_1[MTPL_Ind_i_t, ])))
# Mean [median] risk premium
rbind(c(sum(MTPL_N_X_Sev[, 2]), MTPL_Exp_Freq %*% (lamb_1 * mu_phi_1), 
        MTPL_Exp_Freq %*% (lamb_post_1 * mu_phi_post_1)) / sum(MTPL_Exp_Freq),
      c(weightedMedian(x = Sizes_N / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = lamb_1 * mu_phi_1, w = MTPL_Exp_Freq), 
        weightedMedian(x = lamb_post_1 * mu_phi_post_1, w = MTPL_Exp_Freq)))

##### Predictions K=2 #####
### Estimated prior parameters
cbind(tail(Optim_2$theta_N, 1), tail(Optim_2$theta_X, 1))

### Estimated prior transition probabilities
100 * Optim_2$theta_Z

### Mean predictors
# Prior
lamb_2 <- pmin(pmax(MTPL_Exp_Freq * exp(MTPL_RF_Freq %*% Optim_2$theta_N[
  1:ncol(MTPL_RF_Freq), ]), LB_eps), 1e+15) / MTPL_Exp_Freq
mu_phi_2 <- t(pmin(pmax(Optim_2$theta_X[1, ] * t(exp(MTPL_RF_Sev_Freq %*% Optim_2$theta_X[
  2:(1 + ncol(MTPL_RF_Sev)), ])), LB_eps), 1e+15) / Optim_2$theta_X[1, ])
# Posterior
expo_lamb_2 <- pmin(pmax(MTPL_Exp_Freq * lamb_2, LB_eps), 1e+15)
mu_2 <- pmin(pmax(t(Optim_2$theta_X[1, ] * t(mu_phi_2)), LB_eps), 1e+15)
Update_alpha_U_2 <- rep(0, nrow(MTPL_RF_Freq))
Update_beta_U_2 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_alpha_V_2 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_beta_V_2 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_2 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_2[unique(MTPL_Ind_i_t), ] <- as.matrix(aggregate(
  mu_2[MTPL_Ind_i_t, ] ~ MTPL_Ind_i_t, FUN = sum)[, -1])
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
for (t in 2:max(MTPL_Ind_t_Ti[, 2])){
  t_Current <- which(MTPL_Ind_t_Ti[, 1] == t)
  t_Previous <- which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t))
  Update_alpha_U_2[t_Current] <- Update_alpha_U_2[t_Previous] + MTPL_N_Freq[t_Previous]
  Update_beta_U_2[t_Current, ] <- Update_beta_U_2[t_Previous, ] + expo_lamb_2[t_Previous, ]
  Update_alpha_V_2[t_Current, ] <- Update_alpha_V_2[t_Previous, ] + alpha_N_2[t_Previous, ]
  Update_beta_V_2[t_Current, ] <- Update_beta_V_2[t_Previous, ] + Sizes_N[t_Previous] %*% 
    t(Optim_2$theta_X[1, ])
}
lamb_post_2 <- lamb_2 * sapply(1:K, function(j) 
  (tail(Optim_2$theta_N, 1)[j] + Update_alpha_U_2) /
    (tail(Optim_2$theta_N, 1)[j] + Update_beta_U_2[, j]))
mu_phi_post_2 <- mu_phi_2 * sapply(1:K, function(j) 
  (tail(Optim_2$theta_X, 1)[j] - 1 + Update_beta_V_2[, j]) /
    (tail(Optim_2$theta_X, 1)[j] - 1 + Update_alpha_V_2[, j]))

### Profile assignment probabilities
# Prior
Prob_Z_prior_2 <- matrix(NA, nrow = nrow(MTPL_RF_Freq), ncol = K)
Prob_Z_prior_2[which(MTPL_Ind_t_Ti[, 1] == 1), ] <- Optim_2$theta_Z[1, ]
for (t in 2:max(MTPL_Ind_t_Ti[, 2])) {
  Prob_Z_prior_2[which(MTPL_Ind_t_Ti[, 1] == t), ] <- Prob_Z_prior_2[
    which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t)), ] %*%
    Optim_2$theta_Z[-1, ]
}
# Posterior
Prob_Z_post_2 <- Optim_2$Prob_Z_post_RP

### Response predictions
# Mean claim frequency
rbind(c(sum(MTPL_N_Freq), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_2 * lamb_2, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_2 * lamb_post_2, 1, sum)) / sum(MTPL_Exp_Freq),
      c(weightedMedian(x = MTPL_N_Freq / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_2 * lamb_2, 1, sum), w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_post_2 * lamb_post_2, 1, sum), w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:2, function(j) 
        c(MTPL_Exp_Freq %*% lamb_2[, j] / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_2[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(1:2, function(j) 
            c(MTPL_Exp_Freq %*% lamb_post_2[, j] / sum(MTPL_Exp_Freq), weightedMedian(
              x = lamb_post_2[, j], w = MTPL_Exp_Freq))))))
# Mean [median] claim severity
rbind(c(mean(MTPL_X_Sev), 
        mean(apply(Prob_Z_prior_2[MTPL_Ind_i_t, ] * mu_phi_2[MTPL_Ind_i_t, ], 1, sum)),
        mean(apply(Prob_Z_post_2[MTPL_Ind_i_t, ] * mu_phi_post_2[MTPL_Ind_i_t, ], 1, sum))),
      c(median(MTPL_X_Sev),
        median(apply(Prob_Z_prior_2[MTPL_Ind_i_t, ] * mu_phi_2[MTPL_Ind_i_t, ], 1, sum)),
        median(apply(Prob_Z_post_2[MTPL_Ind_i_t, ] * mu_phi_post_2[MTPL_Ind_i_t, ], 1, sum))),
      cbind(NA, as.vector(sapply(1:2, function(j) 
        c(mean(mu_phi_2[MTPL_Ind_i_t, j]), median(mu_phi_2[MTPL_Ind_i_t, j])))), 
        as.vector(sapply(1:2, function(j) 
          c(mean(mu_phi_post_2[MTPL_Ind_i_t, j]), median(mu_phi_post_2[MTPL_Ind_i_t, j]))))))
# Mean [median] risk premium
rbind(c(sum(MTPL_N_X_Sev[, 2]), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_2 * lamb_2 * mu_phi_2, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_2 * lamb_post_2 * mu_phi_post_2, 1, sum)) / 
        sum(MTPL_Exp_Freq),
      c(weightedMedian(x = Sizes_N / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_2 * lamb_2 * mu_phi_2, 1, sum), 
                       w = MTPL_Exp_Freq), weightedMedian(x = apply(
                         Prob_Z_post_2 * lamb_post_2 * mu_phi_post_2, 1, sum), 
                         w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:2, function(j) 
        c(MTPL_Exp_Freq %*% (lamb_2[, j] * mu_phi_2[, j]) / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_2[, j] * mu_phi_2[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(
            1:2, function(j) c(MTPL_Exp_Freq %*% (lamb_post_2[, j] * mu_phi_post_2[, j]) / 
                                 sum(MTPL_Exp_Freq), weightedMedian(
              x = lamb_post_2[, j] * mu_phi_post_2[, j], w = MTPL_Exp_Freq))))))

### Estimated covariances
# Prior
summary(apply(Prob_Z_prior_2 * lamb_2 * (mu_phi_2 - apply(
  Prob_Z_prior_2 * mu_phi_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_prior_2 * mu_phi_2 * (lamb_2 - apply(
  Prob_Z_prior_2 * lamb_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_prior_2 * lamb_2 * (mu_phi_2 - apply(
  Prob_Z_prior_2 * mu_phi_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))
# Posterior
summary(apply(Prob_Z_post_2 * lamb_post_2 * (mu_phi_post_2 - apply(
  Prob_Z_post_2 * mu_phi_post_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_post_2 * mu_phi_post_2 * (lamb_post_2 - apply(
  Prob_Z_post_2 * lamb_post_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_post_2 * lamb_post_2 * (mu_phi_post_2 - apply(
  Prob_Z_post_2 * mu_phi_post_2, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))

##### Predictions K=3 #####
### Estimated prior parameters
cbind(tail(Optim_3$theta_N, 1), tail(Optim_3$theta_X, 1))

### Estimated prior transition probabilities
100 * Optim_3$theta_Z

### Mean predictors
# Prior
lamb_3 <- pmin(pmax(MTPL_Exp_Freq * exp(MTPL_RF_Freq %*% Optim_3$theta_N[
  1:ncol(MTPL_RF_Freq), ]), LB_eps), 1e+15) / MTPL_Exp_Freq
mu_phi_3 <- t(pmin(pmax(Optim_3$theta_X[1, ] * t(exp(MTPL_RF_Sev_Freq %*% Optim_3$theta_X[
  2:(1 + ncol(MTPL_RF_Sev)), ])), LB_eps), 1e+15) / Optim_3$theta_X[1, ])
# Posterior
expo_lamb_3 <- pmin(pmax(MTPL_Exp_Freq * lamb_3, LB_eps), 1e+15)
mu_3 <- pmin(pmax(t(Optim_3$theta_X[1, ] * t(mu_phi_3)), LB_eps), 1e+15)
Update_alpha_U_3 <- rep(0, nrow(MTPL_RF_Freq))
Update_beta_U_3 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_alpha_V_3 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_beta_V_3 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_3 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_3[unique(MTPL_Ind_i_t), ] <- as.matrix(aggregate(
  mu_3[MTPL_Ind_i_t, ] ~ MTPL_Ind_i_t, FUN = sum)[, -1])
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
for (t in 2:max(MTPL_Ind_t_Ti[, 2])){
  t_Current <- which(MTPL_Ind_t_Ti[, 1] == t)
  t_Previous <- which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t))
  Update_alpha_U_3[t_Current] <- Update_alpha_U_3[t_Previous] + MTPL_N_Freq[t_Previous]
  Update_beta_U_3[t_Current, ] <- Update_beta_U_3[t_Previous, ] + expo_lamb_3[t_Previous, ]
  Update_alpha_V_3[t_Current, ] <- Update_alpha_V_3[t_Previous, ] + alpha_N_3[t_Previous, ]
  Update_beta_V_3[t_Current, ] <- Update_beta_V_3[t_Previous, ] + Sizes_N[t_Previous] %*% 
    t(Optim_3$theta_X[1, ])
}
lamb_post_3 <- lamb_3 * sapply(1:K, function(j) 
  (tail(Optim_3$theta_N, 1)[j] + Update_alpha_U_3) /
    (tail(Optim_3$theta_N, 1)[j] + Update_beta_U_3[, j]))
mu_phi_post_3 <- mu_phi_3 * sapply(1:K, function(j) 
  (tail(Optim_3$theta_X, 1)[j] - 1 + Update_beta_V_3[, j]) /
    (tail(Optim_3$theta_X, 1)[j] - 1 + Update_alpha_V_3[, j]))

### Profile assignment probabilities
# Prior
Prob_Z_prior_3 <- matrix(NA, nrow = nrow(MTPL_RF_Freq), ncol = K)
Prob_Z_prior_3[which(MTPL_Ind_t_Ti[, 1] == 1), ] <- Optim_3$theta_Z[1, ]
for (t in 2:max(MTPL_Ind_t_Ti[, 2])) {
  Prob_Z_prior_3[which(MTPL_Ind_t_Ti[, 1] == t), ] <- Prob_Z_prior_3[
    which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t)), ] %*%
    Optim_3$theta_Z[-1, ]
}
# Posterior
Prob_Z_post_3 <- Optim_3$Prob_Z_post_RP

### Response predictions
# Mean claim frequency
rbind(c(sum(MTPL_N_Freq), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_3 * lamb_3, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_3 * lamb_post_3, 1, sum)) / sum(MTPL_Exp_Freq),
      c(weightedMedian(x = MTPL_N_Freq / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_3 * lamb_3, 1, sum), w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_post_3 * lamb_post_3, 1, sum), w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:3, function(j) 
        c(MTPL_Exp_Freq %*% lamb_3[, j] / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_3[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(1:3, function(j) 
            c(MTPL_Exp_Freq %*% lamb_post_3[, j] / sum(MTPL_Exp_Freq), weightedMedian(
              x = lamb_post_3[, j], w = MTPL_Exp_Freq))))))
# Mean [median] claim severity
rbind(c(mean(MTPL_X_Sev), 
        mean(apply(Prob_Z_prior_3[MTPL_Ind_i_t, ] * mu_phi_3[MTPL_Ind_i_t, ], 1, sum)),
        mean(apply(Prob_Z_post_3[MTPL_Ind_i_t, ] * mu_phi_post_3[MTPL_Ind_i_t, ], 1, sum))),
      c(median(MTPL_X_Sev),
        median(apply(Prob_Z_prior_3[MTPL_Ind_i_t, ] * mu_phi_3[MTPL_Ind_i_t, ], 1, sum)),
        median(apply(Prob_Z_post_3[MTPL_Ind_i_t, ] * mu_phi_post_3[MTPL_Ind_i_t, ], 1, sum))),
      cbind(NA, as.vector(sapply(1:3, function(j) 
        c(mean(mu_phi_3[MTPL_Ind_i_t, j]), median(mu_phi_3[MTPL_Ind_i_t, j])))), 
        as.vector(sapply(1:3, function(j) 
          c(mean(mu_phi_post_3[MTPL_Ind_i_t, j]), median(mu_phi_post_3[MTPL_Ind_i_t, j]))))))
# Mean [median] risk premium
rbind(c(sum(MTPL_N_X_Sev[, 2]), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_3 * lamb_3 * mu_phi_3, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_3 * lamb_post_3 * mu_phi_post_3, 1, sum)) / 
        sum(MTPL_Exp_Freq),
      c(weightedMedian(x = Sizes_N / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_3 * lamb_3 * mu_phi_3, 1, sum), 
                       w = MTPL_Exp_Freq), weightedMedian(x = apply(
                         Prob_Z_post_3 * lamb_post_3 * mu_phi_post_3, 1, sum), 
                         w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:3, function(j) 
        c(MTPL_Exp_Freq %*% (lamb_3[, j] * mu_phi_3[, j]) / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_3[, j] * mu_phi_3[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(
            1:3, function(j) c(MTPL_Exp_Freq %*% (lamb_post_3[, j] * mu_phi_post_3[, j]) / 
                                 sum(MTPL_Exp_Freq), weightedMedian(
                                   x = lamb_post_3[, j] * mu_phi_post_3[, j], w = MTPL_Exp_Freq))))))

### Estimated covariances
# Prior
summary(apply(Prob_Z_prior_3 * lamb_3 * (mu_phi_3 - apply(
  Prob_Z_prior_3 * mu_phi_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_prior_3 * mu_phi_3 * (lamb_3 - apply(
  Prob_Z_prior_3 * lamb_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_prior_3 * lamb_3 * (mu_phi_3 - apply(
  Prob_Z_prior_3 * mu_phi_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))
# Posterior
summary(apply(Prob_Z_post_3 * lamb_post_3 * (mu_phi_post_3 - apply(
  Prob_Z_post_3 * mu_phi_post_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_post_3 * mu_phi_post_3 * (lamb_post_3 - apply(
  Prob_Z_post_3 * lamb_post_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_post_3 * lamb_post_3 * (mu_phi_post_3 - apply(
  Prob_Z_post_3 * mu_phi_post_3, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))

##### Predictions K=4 #####
### Estimated prior parameters
cbind(tail(Optim_4$theta_N, 1), tail(Optim_4$theta_X, 1))

### Estimated prior transition probabilities
100 * Optim_4$theta_Z

### Mean predictors
# Prior
lamb_4 <- pmin(pmax(MTPL_Exp_Freq * exp(MTPL_RF_Freq %*% Optim_4$theta_N[
  1:ncol(MTPL_RF_Freq), ]), LB_eps), 1e+15) / MTPL_Exp_Freq
mu_phi_4 <- t(pmin(pmax(Optim_4$theta_X[1, ] * t(exp(MTPL_RF_Sev_Freq %*% Optim_4$theta_X[
  2:(1 + ncol(MTPL_RF_Sev)), ])), LB_eps), 1e+15) / Optim_4$theta_X[1, ])
# Posterior
expo_lamb_4 <- pmin(pmax(MTPL_Exp_Freq * lamb_4, LB_eps), 1e+15)
mu_4 <- pmin(pmax(t(Optim_4$theta_X[1, ] * t(mu_phi_4)), LB_eps), 1e+15)
Update_alpha_U_4 <- rep(0, nrow(MTPL_RF_Freq))
Update_beta_U_4 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_alpha_V_4 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
Update_beta_V_4 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_4 <- matrix(0, nrow = nrow(MTPL_RF_Freq), ncol = K)
alpha_N_4[unique(MTPL_Ind_i_t), ] <- as.matrix(aggregate(
  mu_4[MTPL_Ind_i_t, ] ~ MTPL_Ind_i_t, FUN = sum)[, -1])
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
for (t in 2:max(MTPL_Ind_t_Ti[, 2])){
  t_Current <- which(MTPL_Ind_t_Ti[, 1] == t)
  t_Previous <- which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t))
  Update_alpha_U_4[t_Current] <- Update_alpha_U_4[t_Previous] + MTPL_N_Freq[t_Previous]
  Update_beta_U_4[t_Current, ] <- Update_beta_U_4[t_Previous, ] + expo_lamb_4[t_Previous, ]
  Update_alpha_V_4[t_Current, ] <- Update_alpha_V_4[t_Previous, ] + alpha_N_4[t_Previous, ]
  Update_beta_V_4[t_Current, ] <- Update_beta_V_4[t_Previous, ] + Sizes_N[t_Previous] %*% 
    t(Optim_4$theta_X[1, ])
}
lamb_post_4 <- lamb_4 * sapply(1:K, function(j) 
  (tail(Optim_4$theta_N, 1)[j] + Update_alpha_U_4) /
    (tail(Optim_4$theta_N, 1)[j] + Update_beta_U_4[, j]))
mu_phi_post_4 <- mu_phi_4 * sapply(1:K, function(j) 
  (tail(Optim_4$theta_X, 1)[j] - 1 + Update_beta_V_4[, j]) /
    (tail(Optim_4$theta_X, 1)[j] - 1 + Update_alpha_V_4[, j]))

### Profile assignment probabilities
# Prior
Prob_Z_prior_4 <- matrix(NA, nrow = nrow(MTPL_RF_Freq), ncol = K)
Prob_Z_prior_4[which(MTPL_Ind_t_Ti[, 1] == 1), ] <- Optim_4$theta_Z[1, ]
for (t in 2:max(MTPL_Ind_t_Ti[, 2])) {
  Prob_Z_prior_4[which(MTPL_Ind_t_Ti[, 1] == t), ] <- Prob_Z_prior_4[
    which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t)), ] %*%
    Optim_4$theta_Z[-1, ]
}
# Posterior
Prob_Z_post_4 <- Optim_4$Prob_Z_post_RP

### Response predictions
# Mean claim frequency
rbind(c(sum(MTPL_N_Freq), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_4 * lamb_4, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_4 * lamb_post_4, 1, sum)) / sum(MTPL_Exp_Freq),
      c(weightedMedian(x = MTPL_N_Freq / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_4 * lamb_4, 1, sum), w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_post_4 * lamb_post_4, 1, sum), w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:4, function(j) 
        c(MTPL_Exp_Freq %*% lamb_4[, j] / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_4[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(1:4, function(j) 
            c(MTPL_Exp_Freq %*% lamb_post_4[, j] / sum(MTPL_Exp_Freq), weightedMedian(
              x = lamb_post_4[, j], w = MTPL_Exp_Freq))))))
# Mean [median] claim severity
rbind(c(mean(MTPL_X_Sev), 
        mean(apply(Prob_Z_prior_4[MTPL_Ind_i_t, ] * mu_phi_4[MTPL_Ind_i_t, ], 1, sum)),
        mean(apply(Prob_Z_post_4[MTPL_Ind_i_t, ] * mu_phi_post_4[MTPL_Ind_i_t, ], 1, sum))),
      c(median(MTPL_X_Sev),
        median(apply(Prob_Z_prior_4[MTPL_Ind_i_t, ] * mu_phi_4[MTPL_Ind_i_t, ], 1, sum)),
        median(apply(Prob_Z_post_4[MTPL_Ind_i_t, ] * mu_phi_post_4[MTPL_Ind_i_t, ], 1, sum))),
      cbind(NA, as.vector(sapply(1:4, function(j) 
        c(mean(mu_phi_4[MTPL_Ind_i_t, j]), median(mu_phi_4[MTPL_Ind_i_t, j])))), 
        as.vector(sapply(1:4, function(j) 
          c(mean(mu_phi_post_4[MTPL_Ind_i_t, j]), median(mu_phi_post_4[MTPL_Ind_i_t, j]))))))
# Mean [median] risk premium
rbind(c(sum(MTPL_N_X_Sev[, 2]), 
        MTPL_Exp_Freq %*% apply(Prob_Z_prior_4 * lamb_4 * mu_phi_4, 1, sum),
        MTPL_Exp_Freq %*% apply(Prob_Z_post_4 * lamb_post_4 * mu_phi_post_4, 1, sum)) / 
        sum(MTPL_Exp_Freq),
      c(weightedMedian(x = Sizes_N / MTPL_Exp_Freq, w = MTPL_Exp_Freq),
        weightedMedian(x = apply(Prob_Z_prior_4 * lamb_4 * mu_phi_4, 1, sum), 
                       w = MTPL_Exp_Freq), weightedMedian(x = apply(
                         Prob_Z_post_4 * lamb_post_4 * mu_phi_post_4, 1, sum), 
                         w = MTPL_Exp_Freq)),
      cbind(NA, as.vector(sapply(1:4, function(j) 
        c(MTPL_Exp_Freq %*% (lamb_4[, j] * mu_phi_4[, j]) / sum(MTPL_Exp_Freq), weightedMedian(
          x = lamb_4[, j] * mu_phi_4[, j], w = MTPL_Exp_Freq)))), as.vector(sapply(
            1:4, function(j) c(MTPL_Exp_Freq %*% (lamb_post_4[, j] * mu_phi_post_4[, j]) / 
                                 sum(MTPL_Exp_Freq), weightedMedian(
                                   x = lamb_post_4[, j] * mu_phi_post_4[, j], w = MTPL_Exp_Freq))))))

### Estimated covariances
# Prior
summary(apply(Prob_Z_prior_4 * lamb_4 * (mu_phi_4 - apply(
  Prob_Z_prior_4 * mu_phi_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_prior_4 * mu_phi_4 * (lamb_4 - apply(
  Prob_Z_prior_4 * lamb_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_prior_4 * lamb_4 * (mu_phi_4 - apply(
  Prob_Z_prior_4 * mu_phi_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))
# Posterior
summary(apply(Prob_Z_post_4 * lamb_post_4 * (mu_phi_post_4 - apply(
  Prob_Z_post_4 * mu_phi_post_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
summary(apply(Prob_Z_post_4 * mu_phi_post_4 * (lamb_post_4 - apply(
  Prob_Z_post_4 * lamb_post_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)])
c(sum(MTPL_N_Freq > 0), sum(apply(Prob_Z_post_4 * lamb_post_4 * (mu_phi_post_4 - apply(
  Prob_Z_post_4 * mu_phi_post_4, 1, sum)), 1, sum)[which(MTPL_N_Freq > 0)] > 0))

##### Performance #####
### Log-likelihood
logL_1 <- tail(Optim_1$ell_BW, 1)
logL_2 <- tail(Optim_2$ell_BW, 1)
logL_3 <- tail(Optim_3$ell_BW, 1)
logL_4 <- tail(Optim_4$ell_BW, 1)

### Number of parameters
P_1 <- length(Optim_1$theta_N) + length(Optim_1$theta_X)
P_2 <- length(Optim_2$theta_Z) + length(Optim_2$theta_N) + length(Optim_2$theta_X)
P_3 <- length(Optim_3$theta_Z) + length(Optim_3$theta_N) + length(Optim_3$theta_X)
P_4 <- length(Optim_4$theta_Z) + length(Optim_4$theta_N) + length(Optim_4$theta_X)

### Akaike Information Criterion
AIC_1 <- 2 * P_1 - 2 * logL_1
AIC_2 <- 2 * P_2 - 2 * logL_2
AIC_3 <- 2 * P_3 - 2 * logL_3
AIC_4 <- 2 * P_4 - 2 * logL_4

### Bayesian Information Criterion
BIC_1 <- log(nrow(MTPL_RF_Freq)) * P_1 - 2 * logL_1
BIC_2 <- log(nrow(MTPL_RF_Freq)) * P_2 - 2 * logL_2
BIC_3 <- log(nrow(MTPL_RF_Freq)) * P_3 - 2 * logL_3
BIC_4 <- log(nrow(MTPL_RF_Freq)) * P_4 - 2 * logL_4

### Prior loss ratio
LR_prior_1 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * lamb_1 * mu_phi_1)
LR_prior_2 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_prior_2 * lamb_2 * mu_phi_2, 1, sum))
LR_prior_3 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_prior_3 * lamb_3 * mu_phi_3, 1, sum))
LR_prior_4 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_prior_4 * lamb_4 * mu_phi_4, 1, sum))

### Posterior loss ratio
LR_post_1 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * lamb_post_1 * mu_phi_post_1)
LR_post_2 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_post_2 * lamb_post_2 * mu_phi_post_2, 1, sum))
LR_post_3 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_post_3 * lamb_post_3 * mu_phi_post_3, 1, sum))
LR_post_4 <- 100 * sum(MTPL_X_Sev) / sum(MTPL_Exp_Freq * apply(
  Prob_Z_post_4 * lamb_post_4 * mu_phi_post_4, 1, sum))

##### Ordered Lorenz curves #####
### Function for calculating the ordered Lorenz curve
Lorenz_H <- function(H, Ref, Act, Pred) {
  # Order based on the predictions relative to the benchmark
  OLC <- as.matrix(cbind(Ref, Act, Pred / Ref))
  OLC <- OLC[order(OLC[, 3], decreasing = FALSE), ]
  
  # Determine the empirical cumulative distribution function
  OLC <- cbind(cumsum(OLC[, 2]) / sum(OLC[, 2]), cumsum(OLC[, 1]) / sum(OLC[, 1]), OLC)
  Rep <- sapply(1:H, function(z) which((OLC[, 2] > ((z - 1) / H)) & (OLC[, 2] <= (z / H))))
  OLC_H <- as.matrix(cbind(c(0, sapply(1:H, function(z) tail(OLC[Rep[[z]], 1], 1))), 0:H))
  
  # Return the ordered Lorenz curve
  return(OLC_H)
}

### Calculate the ordered Lorenz curves
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
H <- 1000
Lorenz_prior_1 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                             MTPL_Exp_Freq * lamb_1 * mu_phi_1)
Lorenz_post_1 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                            MTPL_Exp_Freq * lamb_post_1 * mu_phi_post_1)
Lorenz_prior_2 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                             MTPL_Exp_Freq * apply(Prob_Z_prior_2 * lamb_2 * mu_phi_2, 1, sum))
Lorenz_post_2 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                            MTPL_Exp_Freq * apply(Prob_Z_post_2 * lamb_post_2 * mu_phi_post_2, 1, sum))
Lorenz_prior_3 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                             MTPL_Exp_Freq * apply(Prob_Z_prior_3 * lamb_3 * mu_phi_3, 1, sum))
Lorenz_post_3 <- Lorenz_H(H = H, Ref = MTPL_Exp_Freq, Act = Sizes_N, Pred = 
                            MTPL_Exp_Freq * apply(Prob_Z_post_3 * lamb_post_3 * mu_phi_post_3, 1, sum))
Lorenz_pp_1_2_3 <- data.frame(rbind(Lorenz_prior_1[, c(1, 2)], Lorenz_post_1[, c(1, 2)], 
                                    Lorenz_prior_2[, c(1, 2)], Lorenz_post_2[, c(1, 2)],
                                    Lorenz_prior_3[, c(1, 2)], Lorenz_post_3[, c(1, 2)]))

##### Ratio Gini coefficients #####
### Load function for calculating ratio Gini coefficient
rGini <- dget('.../Common_functions/rGini.R')

### Construction of model forecasts
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
Mod_fc <- cbind(
  MTPL_Exp_Freq * lamb_1 * alp_bet_1,
  MTPL_Exp_Freq * lamb_post_1 * alp_bet_post_1,
  MTPL_Exp_Freq * apply(Prob_Z_prior_2 * lamb_2 * alp_bet_2, 1, sum),
  MTPL_Exp_Freq * apply(Prob_Z_post_2 * lamb_post_2 * alp_bet_post_2, 1, sum),
  MTPL_Exp_Freq * apply(Prob_Z_prior_3 * lamb_3 * alp_bet_3, 1, sum),
  MTPL_Exp_Freq * apply(Prob_Z_post_3 * lamb_post_3 * alp_bet_post_3, 1, sum),
  MTPL_Exp_Freq * apply(Prob_Z_prior_4 * lamb_4 * alp_bet_4, 1, sum),
  MTPL_Exp_Freq * apply(Prob_Z_post_4 * lamb_post_4 * alp_bet_post_4, 1, sum)
)

### Calculation of ratio Gini coefficients
ratio_Gini_coef <- matrix(NA, nrow = ncol(Mod_fc), ncol = ncol(Mod_fc))
# Loop over all benchmark models:
for (b in 1:ncol(Mod_fc)) {
  print(paste('Benchmark: ', b, '/', ncol(Mod_fc), sep = ''))
  # Loop over all alternative models:
  if (b == 1) {
    a_b <- c(3, 5, 7)
  } else if (b == 2) {
    a_b <- c(4, 6, 8)
  } else if (b %in% c(3, 5, 7)) {
    a_b <- 1
  } else {
    a_b <- 2
  }
  for (a in a_b) {
    print(paste('   Alternative: ', a, '/', ncol(Mod_fc), sep = ''))
    # Determine the ratio Gini coefficient:
    ratio_Gini_coef[b, a] <- 100 * rGini(Candidate = Mod_fc[, a], Reference = Mod_fc[, b], 
                                         Actual = Sizes_N)[[1]][1]
  }
}

##### Spearman's rho & Kendall's tau #####
### Load function for calculating Spearman's rho and Kendall's tau
rho_tau <- dget('.../Common_functions/rho_tau.R')

### K = 2 (prior)
# Input
K <- 2
size <- tail(Optim_2$theta_N, 1)
prob <- t(as.vector(tail(Optim_2$theta_N, 1)) / (
  as.vector(tail(Optim_2$theta_N, 1)) + t(lamb_2)))
shape1 <- 1
shape2 <- mu_2_select
shape3 <- tail(Optim_2$theta_X, 1)
scale <- (tail(Optim_2$theta_X, 1) - 1) / Optim_2$theta_X[1, ]
Prob_Z <- Prob_Z_prior_2
# Output
rho_prior_2 <- c()
tau_prior_2 <- c()
for (i in 1:nrow(Prob_Z)) {
  dep_i <- rho_tau(K = K, size = as.vector(size), prob = as.vector(prob[i, ]), 
                   shape1 = shape1, shape2 = as.vector(shape2[i, ]), shape3 = as.vector(shape3), 
                   scale = as.vector(scale), Prob_Z = as.vector(Prob_Z[i, ]))
  rho_prior_2[i] <- dep_i[1]
  tau_prior_2[i] <- dep_i[2]
}

### K = 2 (posterior)
K <- 2
size <- sapply(1:K, function(j) tail(Optim_2$theta_N, 1)[j] + Update_alpha_U_2)
prob <- t((as.vector(tail(Optim_2$theta_N, 1)) + t(Update_beta_U_2)) / (
  (as.vector(tail(Optim_2$theta_N, 1)) + t(Update_beta_U_2)) + t(lamb_2)))
shape1 <- 1
shape2 <- mu_2
shape3 <- t(as.vector(tail(Optim_2$theta_X, 1)) + t(Update_alpha_V_2))
scale <- t((as.vector(tail(Optim_2$theta_X, 1)) - 1 + t(Update_beta_V_2)) / Optim_2$theta_X[1, ])
Prob_Z <- Prob_Z_post_2
# Output
rho_post_2 <- c()
tau_post_2 <- c()
for (i in 1:nrow(Prob_Z)) {
  dep_i <- rho_tau(K = K, size = as.vector(size[i, ]), prob = as.vector(prob[i, ]), 
                   shape1 = shape1, shape2 = as.vector(shape2[i, ]), shape3 = as.vector(shape3[i, ]), 
                   scale = as.vector(scale[i, ]), Prob_Z = as.vector(Prob_Z[i, ]))
  rho_post_2[i] <- dep_i[1]
  tau_post_2[i] <- dep_i[2]
}

### K = 3 (prior)
K <- 3
size <- tail(Optim_3$theta_N, 1)
prob <- t(as.vector(tail(Optim_3$theta_N, 1)) / (
  as.vector(tail(Optim_3$theta_N, 1)) + t(lamb_3)))
shape1 <- 1
shape2 <- mu_3_select
shape3 <- tail(Optim_3$theta_X, 1)
scale <- (tail(Optim_3$theta_X, 1) - 1) / Optim_3$theta_X[1, ]
Prob_Z <- Prob_Z_prior_3
# Output
rho_prior_3 <- c()
tau_prior_3 <- c()
for (i in 1:nrow(Prob_Z)) {
  dep_i <- rho_tau(K = K, size = as.vector(size), prob = as.vector(prob[i, ]), 
                   shape1 = shape1, shape2 = as.vector(shape2[i, ]), shape3 = as.vector(shape3), 
                   scale = as.vector(scale), Prob_Z = as.vector(Prob_Z[i, ]))
  rho_prior_3[i] <- dep_i[1]
  tau_prior_3[i] <- dep_i[2]
}

### K = 3 (posterior)
K <- 3
size <- sapply(1:K, function(j) tail(Optim_3$theta_N, 1)[j] + Update_alpha_U_3)
prob <- t((as.vector(tail(Optim_3$theta_N, 1)) + t(Update_beta_U_3)) / (
  (as.vector(tail(Optim_3$theta_N, 1)) + t(Update_beta_U_3)) + t(lamb_3)))
shape1 <- 1
shape2 <- mu_3
shape3 <- t(as.vector(tail(Optim_3$theta_X, 1)) + t(Update_alpha_V_3))
scale <- t((as.vector(tail(Optim_3$theta_X, 1)) - 1 + t(Update_beta_V_3)) / Optim_3$theta_X[1, ])
Prob_Z <- Prob_Z_post_3
# Output
rho_post_3 <- c()
tau_post_3 <- c()
for (i in 1:nrow(Prob_Z)) {
  dep_i <- rho_tau(K = K, size = as.vector(size[i, ]), prob = as.vector(prob[i, ]), 
                   shape1 = shape1, shape2 = as.vector(shape2[i, ]), shape3 = as.vector(shape3[i, ]), 
                   scale = as.vector(scale[i, ]), Prob_Z = as.vector(Prob_Z[i, ]))
  rho_post_3[i] <- dep_i[1]
  tau_post_3[i] <- dep_i[2]
}

##### Underlying risk profiles #####
### Aggregate past claims experience
MTPL_Hist_N_X <- matrix(0, nrow = length(MTPL_N_Freq), ncol = 4)
Sizes_N <- rep(0, nrow(MTPL_RF_Freq))
Sizes_N[unique(MTPL_Ind_i_t)] <- as.vector(aggregate(
  MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)[, -1])
for (t in 2:max(MTPL_Ind_t_Ti[, 2])){
  t_Current <- which(MTPL_Ind_t_Ti[, 1] == t)
  t_Previous <- which((MTPL_Ind_t_Ti[, 1] == (t - 1)) & (MTPL_Ind_t_Ti[, 2] >= t))
  MTPL_Hist_N_X[t_Current, 1] <- MTPL_Hist_N_X[t_Previous, 1] + MTPL_Exp_Freq[t_Previous]
  MTPL_Hist_N_X[t_Current, 2] <- MTPL_Hist_N_X[t_Previous, 2] + MTPL_N_Freq[t_Previous]
  MTPL_Hist_N_X[t_Current, 3] <- MTPL_Hist_N_X[t_Previous, 3] + Sizes_N[t_Previous]
}
MTPL_Hist_N_X[, 4] <- Sizes_N
MTPL_Hist_N_X_df <- as.data.frame(MTPL_Hist_N_X)
colnames(MTPL_Hist_N_X_df) <- c('Exposure', 'Count', 'Size', 'Current')
# Filter on all observations with at least some history (exposure > 0)
MTPL_Hist_N_X_df_adj <- MTPL_Hist_N_X_df
for (ex in 0:ceiling(max(MTPL_Hist_N_X_df$Exposure))) {
  MTPL_Hist_N_X_df_adj[which((MTPL_Hist_N_X_df$Exposure > (ex - 1)) & 
                               (MTPL_Hist_N_X_df$Exposure <= ex)), 1] <- ex
}
MTPL_Hist_N_X_df_adj <- MTPL_Hist_N_X_df_adj[which(MTPL_Hist_N_X_df_adj$Exposure > 0), ]

### Distribution of posterior profile assignment probabilities for K=2
MTPL_Hist_N_X_K2 <- MTPL_Hist_N_X_df
MTPL_Hist_N_X_K2$Post_1 <- Prob_Z_post_2[, 1]
MTPL_Hist_N_X_K2$Post_2 <- Prob_Z_post_2[, 2]
# Filter on all observations with at least some history (exposure > 0)
MTPL_Hist_N_X_K2_adj <- MTPL_Hist_N_X_K2
for (ex in 0:ceiling(max(MTPL_Hist_N_X_K2$Exposure))) {
  MTPL_Hist_N_X_K2_adj[which((MTPL_Hist_N_X_K2$Exposure > (ex - 1)) & 
                              (MTPL_Hist_N_X_K2$Exposure <= ex)), 1] <- ex
}
MTPL_Hist_N_X_K2_adj <- MTPL_Hist_N_X_K2_adj[which(MTPL_Hist_N_X_K2_adj$Exposure > 0), ]

### Distribution of posterior profile assignment probabilities for K=2
MTPL_Hist_N_X_K3 <- MTPL_Hist_N_X_df
MTPL_Hist_N_X_K3$Post_1 <- Prob_Z_post_3[, 1]
MTPL_Hist_N_X_K3$Post_2 <- Prob_Z_post_3[, 2]
MTPL_Hist_N_X_K3$Post_3 <- Prob_Z_post_3[, 3]
# Filter on all observations with at least some history (exposure > 0)
MTPL_Hist_N_X_K3_adj <- MTPL_Hist_N_X_K3
for (ex in 0:ceiling(max(MTPL_Hist_N_X_K3$Exposure))) {
  MTPL_Hist_N_X_K3_adj[which((MTPL_Hist_N_X_K3$Exposure > (ex - 1)) & 
                              (MTPL_Hist_N_X_K3$Exposure <= ex)), 1] <- ex
}
MTPL_Hist_N_X_K3_adj <- MTPL_Hist_N_X_K3_adj[which(MTPL_Hist_N_X_K3_adj$Exposure > 0), ]

### Posterior, Bonus-Malus corrections of prior risk premium for K=1
MTPL_BM_1 <- MTPL_Hist_N_X_df
MTPL_BM_1$RP_prior <- lamb_1 * mu_phi_1
MTPL_BM_1$RP_post <- lamb_post_1 * mu_phi_post_1
MTPL_BM_1$RP_corr <- MTPL_BM_1$RP_post / MTPL_BM_1$RP_prior
# Filter on all observations with at least some history (past exposure > 0)
MTPL_BM_1_adj <- MTPL_BM_1
for (ex in 0:ceiling(max(MTPL_BM_1$Exposure))) {
  MTPL_BM_1_adj[which((MTPL_BM_1$Exposure > (ex - 1)) & 
                        (MTPL_BM_1$Exposure <= ex)), 1] <- ex
}
MTPL_BM_1_adj <- MTPL_BM_1_adj[which(MTPL_BM_1_adj$Exposure > 0), ]

### Posterior, Bonus-Malus corrections of prior risk premium for K=2
MTPL_BM_2 <- MTPL_Hist_N_X_df
MTPL_BM_2$RP_prior <- apply(Prob_Z_prior_2 * lamb_2 * mu_phi_2, 1, sum)
MTPL_BM_2$RP_post <- apply(Prob_Z_post_2 * lamb_post_2 * mu_phi_post_2, 1, sum)
MTPL_BM_2$RP_corr <- MTPL_BM_2$RP_post / MTPL_BM_2$RP_prior
# Filter on all observations with at least some history (past exposure > 0)
MTPL_BM_2_adj <- MTPL_BM_2
for (ex in 0:ceiling(max(MTPL_BM_2$Exposure))) {
  MTPL_BM_2_adj[which((MTPL_BM_2$Exposure > (ex - 1)) & 
                        (MTPL_BM_2$Exposure <= ex)), 1] <- ex
}
MTPL_BM_2_adj <- MTPL_BM_2_adj[which(MTPL_BM_2_adj$Exposure > 0), ]

### Posterior, Bonus-Malus corrections of prior risk premium for K=3
MTPL_BM_3 <- MTPL_Hist_N_X_df
MTPL_BM_3$RP_prior <- apply(Prob_Z_prior_3 * lamb_3 * mu_phi_3, 1, sum)
MTPL_BM_3$RP_post <- apply(Prob_Z_post_3 * lamb_post_3 * mu_phi_post_3, 1, sum)
MTPL_BM_3$RP_corr <- MTPL_BM_3$RP_post / MTPL_BM_3$RP_prior
# Filter on all observations with at least some history (past exposure > 0)
MTPL_BM_3_adj <- MTPL_BM_3
for (ex in 0:ceiling(max(MTPL_BM_3$Exposure))) {
  MTPL_BM_3_adj[which((MTPL_BM_3$Exposure > (ex - 1)) & 
                        (MTPL_BM_3$Exposure <= ex)), 1] <- ex
}
MTPL_BM_3_adj <- MTPL_BM_3_adj[which(MTPL_BM_3_adj$Exposure > 0), ]
