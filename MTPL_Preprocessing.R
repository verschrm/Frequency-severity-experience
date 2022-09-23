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
Packages <- c('dplyr', 'MASS', 'gamlss')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load confidential insurance data on Motor Third Party Liability (MTPL) insurance
#   - MTPL_Freq containing information on the claim counts;
#   - MTPL_Sev containing information of the individual claim severities.

##### Preprocessing #####
### Convert data into dataframe
MTPL_Freq_df <- data.frame(MTPL_Freq)
MTPL_Sev_df <- data.frame(MTPL_Sev)

### Order policies and claims of same customer from oldest to newest
MTPL_Freq_df <- MTPL_Freq_df %>% group_by(Cust_ID) %>% arrange(Start_Date, .by_group = TRUE)
MTPL_Sev_df <- MTPL_Sev_df %>% group_by(Cust_ID) %>% arrange(Start_Date, Sev_Date, .by_group = TRUE)

### Extract key variables
MTPL_N_Freq <- MTPL_Freq_df$Count
MTPL_Exp_Freq <- MTPL_Freq_df$Exposure
MTPL_X_Sev <- MTPL_Sev_df$Size

### Extract risk factors, where dummies have been one-hot encoded
MTPL_RF_Freq <- model.matrix(~ 1 + Cust_Age + Cust_Residence + Veh_Age + Veh_BodyDoors + Veh_CatValue +
                               Veh_FuelType + Veh_Mileage + Veh_PowerWeight + Veh_Region + Veh_Weight,
                             MTPL_Freq_df)
MTPL_RF_Sev <- model.matrix(~ 1 + Cust_Age + Cust_Residence + Veh_Age + Veh_BodyDoors + Veh_CatValue +
                              Veh_FuelType + Veh_Mileage + Veh_PowerWeight + Veh_Region + Veh_Weight,
                            MTPL_Sev_df)
MTPL_RF_Sev_Freq <- model.matrix(~ 1 + Cust_Age + Cust_Residence + Veh_Age + Veh_BodyDoors + Veh_CatValue +
                                   Veh_FuelType + Veh_Mileage + Veh_PowerWeight + Veh_Region + Veh_Weight,
                                 MTPL_Freq_df)

### Deduce the observation period and maximum number of periods observed for each customer
MTPL_Ind_t_Ti <- MTPL_Freq_df %>% group_by(Cust_ID) %>% 
  mutate(Ind_t = order(Start_Date), Ind_Ti = length(Count), .keep = 'used')
MTPL_Ind_t_Ti <- MTPL_Ind_t_Ti[, -c(1:3)]

### Match selected severities to frequencies
MTPL_Freq_nonzero <- MTPL_Freq_df
MTPL_Freq_nonzero$Row_num <- as.integer(1:nrow(MTPL_Freq_df))
MTPL_Freq_nonzero <- MTPL_Freq_nonzero[which(MTPL_Freq_nonzero$Count > 0), 
                                       c('Cust_ID', 'Start_Date', 'Row_num')]
MTPL_Sev_df <- MTPL_Sev_df %>%
  group_by(Cust_ID, Start_Date) %>%
  mutate(Freq_row_Select = MTPL_Freq_nonzero$Row_num[which((MTPL_Freq_nonzero$Cust_ID == Cust_ID[1]) & 
                                                            (MTPL_Freq_nonzero$Start_Date == Start_Date[1]))], 
         Sequence_Select = order(Sev_Date), .keep = 'all')
rm(MTPL_Freq_nonzero)
MTPL_Ind_i_t <- MTPL_Sev_df$Freq_row_Select
MTPL_N_X_Sev <- aggregate(MTPL_X_Sev ~ MTPL_Ind_i_t, FUN = sum)

##### GB2 robustness to outliers #####
### Determine interquartile range of claim sizes
MTPL_Q1 <- quantile(MTPL_X_Sev, probs = 0.25)
MTPL_Q3 <- quantile(MTPL_X_Sev, probs = 0.75)
MTPL_IQR <- MTPL_Q3 - MTPL_Q1

### GB2 distribution based on all observed claim sizes
c(length(TPL_X_Sev), mean(TPL_X_Sev), median(TPL_X_Sev))
G_temp <- glm(formula = MTPL_X_Sev ~ 1 + MTPL_RF_Sev, family = Gamma(link = 'log'))
G_temp$shape <- gamma.shape(G_temp)
GB2_temp <- gamlss(formula = MTPL_X_Sev ~ 1, sigma.formula = MTPL_X_Sev ~ 1,
                   nu.formula = MTPL_X_Sev ~ 1 + MTPL_RF_Sev, tau.formula = MTPL_X_Sev ~ 1,
                   family = GB2(mu.link = 'log', sigma.link = 'log',
                                nu.link = 'log', tau.link = 'log'),
                   mu.start = (2 - 1) / G_temp$shape$alpha, sigma.start = 1,
                   nu.start = G_temp$shape$alpha * as.numeric(G_temp$fitted.values),
                   tau.start = 2, sigma.fix = TRUE,
                   control = gamlss.control(c.crit = 1e-8, n.cyc = 100000, trace = FALSE), 
                   i.control = glim.control(cc = 1e-8, cyc = 100000, bf.cyc = 100000, bf.tol = 1e-8))
GB2_temp_mu <- as.numeric(exp(GB2_temp$mu.coefficients))
GB2_temp_sigma <- 1
GB2_temp_nu <- as.vector(exp(MTPL_RF_Sev %*% GB2_temp$nu.coefficients))
GB2_temp_tau <- as.vector(exp(GB2_temp$tau.coefficients))
summary(GB2_temp_nu * GB2_temp_mu / (GB2_temp_tau - 1))

### GB2 distribution of observed claim sizes limited to some maximal claim amount
Limit_out <- c(100000, 50000, 25000, 10000)
for (lim in Limit_out) {
  # GB2 predicted claim sizes after outlier truncation
  GB2_temp <- gamlss(formula = pmin(MTPL_X_Sev, lim) ~ 1, 
                     sigma.formula = pmin(MTPL_X_Sev, lim) ~ 1,
                     nu.formula = pmin(MTPL_X_Sev, lim) ~ -1 + MTPL_RF_Sev, 
                     tau.formula = pmin(MTPL_X_Sev, lim) ~ 1,
                     family = GB2(mu.link = 'log', sigma.link = 'log',
                                  nu.link = 'log', tau.link = 'log'),
                     mu.start = (as.numeric(G_temp$fitted.values) - 1) / G_temp$shape$alpha, 
                     sigma.start = 1, nu.start = G_temp$shape$alpha, tau.start = 2, sigma.fix = TRUE,
                     control = gamlss.control(c.crit = 1e-8, n.cyc = 100000, trace = FALSE), 
                     i.control = glim.control(cc = 1e-8, cyc = 100000, bf.cyc = 100000, bf.tol = 1e-8))
  GB2_temp_mu <- as.numeric(exp(GB2_temp$mu.coefficients))
  GB2_temp_sigma <- 1
  GB2_temp_nu <- as.vector(exp(MTPL_RF_Sev %*% GB2_temp$nu.coefficients))
  GB2_temp_tau <- as.numeric(exp(GB2_temp$tau.coefficients))
  summary(GB2_temp_nu * GB2_temp_mu / (GB2_temp_tau - 1))
  
  # GB2 predicted claim sizes after outlier removal
  GB2_temp <- gamlss(formula = MTPL_X_Sev[which(MTPL_X_Sev <= lim)] ~ 1, 
                     sigma.formula = MTPL_X_Sev[which(MTPL_X_Sev <= lim)] ~ 1,
                     nu.formula = MTPL_X_Sev[which(MTPL_X_Sev <= lim)] ~ -1 + 
                       MTPL_RF_Sev[which(MTPL_X_Sev <= lim), ], 
                     tau.formula = MTPL_X_Sev[which(MTPL_X_Sev <= lim)] ~ 1,
                     family = GB2(mu.link = 'log', sigma.link = 'log',
                                  nu.link = 'log', tau.link = 'log'),
                     mu.start = (as.numeric(G_temp$fitted.values[which(MTPL_X_Sev <= lim)]) - 1) / 
                       G_temp$shape$alpha, sigma.start = 1, nu.start = G_temp$shape$alpha,
                     tau.start = 2, sigma.fix = TRUE,
                     control = gamlss.control(c.crit = 1e-8, n.cyc = 100000, trace = FALSE), 
                     i.control = glim.control(cc = 1e-8, cyc = 100000, bf.cyc = 100000, bf.tol = 1e-8))
  GB2_temp_mu <- as.numeric(exp(GB2_temp$mu.coefficients))
  GB2_temp_sigma <- 1
  GB2_temp_nu <- as.vector(exp(MTPL_RF_Sev[which(MTPL_X_Sev <= lim), ] %*% GB2_temp$nu.coefficients))
  GB2_temp_tau <- as.numeric(exp(GB2_temp$tau.coefficients))
  summary(GB2_temp_nu * GB2_temp_mu / (GB2_temp_tau - 1))
}

### GB2 distribution of observed claim sizes after outlier truncation based on interquartile range
GB2_temp <- gamlss(formula = pmin(pmax(MTPL_Q1 - 1.5 * MTPL_IQR, MTPL_X_Sev), 
                                  MTPL_Q3 + 1.5 * MTPL_IQR) ~ 1, 
                   sigma.formula = pmin(pmax(MTPL_Q1 - 1.5 * MTPL_IQR, MTPL_X_Sev), 
                                        MTPL_Q3 + 1.5 * MTPL_IQR) ~ 1,
                   nu.formula = pmin(pmax(MTPL_Q1 - 1.5 * MTPL_IQR, MTPL_X_Sev), 
                                     MTPL_Q3 + 1.5 * MTPL_IQR) ~ -1 + MTPL_RF_Sev, 
                   tau.formula = pmin(pmax(MTPL_Q1 - 1.5 * MTPL_IQR, MTPL_X_Sev), 
                                      MTPL_Q3 + 1.5 * MTPL_IQR) ~ 1,
                   family = GB2(mu.link = 'log', sigma.link = 'log',
                                nu.link = 'log', tau.link = 'log'),
                   mu.start = (as.numeric(G_temp$fitted.values) - 1) / G_temp$shape$alpha, 
                   sigma.start = 1, nu.start = G_temp$shape$alpha,ntau.start = 2, sigma.fix = TRUE,
                   control = gamlss.control(c.crit = 1e-8, n.cyc = 100000, trace = FALSE), 
                   i.control = glim.control(cc = 1e-8, cyc = 100000, bf.cyc = 100000, bf.tol = 1e-8))
GB2_temp_mu <- as.numeric(exp(GB2_temp$mu.coefficients))
GB2_temp_sigma <- 1
GB2_temp_nu <- as.vector(exp(MTPL_RF_Sev %*% GB2_temp$nu.coefficients))
GB2_temp_tau <- as.numeric(exp(GB2_temp$tau.coefficients))
summary(GB2_temp_nu * GB2_temp_mu / (GB2_temp_tau - 1))

### GB2 distribution of observed claim sizes after outlier removal based on interquartile range
GB2_temp <- gamlss(formula = MTPL_X_Sev[which((MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                                                (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR))] ~ 1, 
                   sigma.formula = MTPL_X_Sev[which((MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                                                      (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR))] ~ 1,
                   nu.formula = MTPL_X_Sev[which((MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                                                   (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR))] ~ -1 + 
                     MTPL_RF_Sev[which((MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                                         (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR)), ], 
                   tau.formula = MTPL_X_Sev[which((MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                                                    (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR))] ~ 1,
                   family = GB2(mu.link = 'log', sigma.link = 'log',
                                 nu.link = 'log', tau.link = 'log'),
                   mu.start = (as.numeric(G_temp$fitted.values[which(
                     (MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & 
                       (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR))]) - 1) / G_temp$shape$alpha, 
                   sigma.start = 1, nu.start = G_temp$shape$alpha, tau.start = 2, sigma.fix = TRUE,
                   control = gamlss.control(c.crit = 1e-8, n.cyc = 100000, trace = FALSE), 
                   i.control = glim.control(cc = 1e-8, cyc = 100000, bf.cyc = 100000, bf.tol = 1e-8))
GB2_temp_mu <- as.numeric(exp(GB2_temp$mu.coefficients))
GB2_temp_sigma <- 1
GB2_temp_nu <- as.vector(exp(MTPL_RF_Sev[which(
  (MTPL_X_Sev >= MTPL_Q1 - 1.5 * MTPL_IQR) & (MTPL_X_Sev <= MTPL_Q3 + 1.5 * MTPL_IQR)), ] %*% 
    GB2_temp$nu.coefficients))
GB2_temp_tau <- as.numeric(exp(GB2_temp$tau.coefficients))
summary(GB2_temp_nu * GB2_temp_mu / (GB2_temp_tau - 1))

rm(G_temp, GB2_temp, GB2_temp_mu, GB2_temp_sigma, GB2_temp_nu, GB2_temp_tau)

### Store the final results
