##
## Libraries and functions
##

Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcMixed.cpp')

library(bayesTreeRing)
library(myFunctions)

##
## load the data
##

load('~/treeRing/datafiles/hudsonValleyData.RData')
species <- as.integer(as.factor(species.crn[ - 1]))
num_species <- length(unique(species))
p <- length(species) 
y <- as.matrix(crn.dat[ - 1, 2:35])
H <- !is.na(y)
y[!H] <- 0L


J <- rep(1, 110)
year_idx <- rep(1:12, 110)

Temp <- t(matrix(base::rowMeans(Temp.avg.dat[ - which(Temp.avg.dat$Year > 2004), 3:36]), 110, 12))
P <- t(log(matrix(base::rowMeans(Precip.dat[ - which(Precip.dat$Year > 2004), 3:36]), 110, 12)))

day <- as.matrix(read.table('~/Mechanistic-Tree-Ring/datafiles/daylength', sep = "", header = TRUE))
day[day == -9999] <- NA
day_len <- apply(day, 2, mean, na.rm = TRUE)
day_len <- day_len / max(day_len)

## fit empirical Bayes climate model
process <- makeMCARTrend(Temp, P)
  
##
## set up priors
##

## calibration model priors
mu_mu_beta_0 <- 1
s2_mu_beta_0 <- 1
alpha_s2_beta_0 <- 1
beta_s2_beta_0 <- 1
mu_mu_beta_1 <- 1
s2_mu_beta_1 <- 1
alpha_s2_beta_1 <- 1
beta_s2_beta_1 <- 1
alpha_s2 <- 1
beta_s2 <- 1

## growth model priors
## Temp_min is beta(9, 5, 0, 9)
alpha_Temp_min <- 9
beta_Temp_min <- 5
Temp_min_lower <- 0
Temp_min_upper <- 9 
## Temp_max is beta(3.5, 3.5, 10, 24)
alpha_Temp_max <- 3.5
beta_Temp_max <- 3.5
Temp_max_lower <- 10
Temp_max_upper <- 24
## P_min is beta(3.5, 3.5, 65, 85)
alpha_P_min <- 3.5
beta_P_min <- 3.5
P_min_lower <- 65
P_min_upper <- 85
## P_max is beta(3.5, 3.5, 85, 105)
alpha_P_max <- 3.5
beta_P_max <- 3.5
P_max_lower <- 85
P_max_upper <- 105
## Probit growth parameters
mu_mu_gamma_T <- 13
s2_mu_gamma_T <- 4
alpha_s2_gamma_T <- 2
beta_s2_gamma_T <- 0.5
mu_mu_xi2_T <- 0
s2_mu_xi2_T <- 1
alpha_s2_xi2_T <- 2
beta_s2_xi2_T <- 0.5
mu_mu_gamma_P <- 85
s2_mu_gamma_P <- 16
alpha_s2_gamma_P <- 2
beta_s2_gamma_P <- 0.5
mu_mu_xi2_P <- 0
s2_mu_xi2_P <- 2
alpha_s2_xi2_P <- 2
beta_s2_xi2_P <- 0.5

## tuning parameters
Temp_min_tune <- 1.0
Temp_max_tune <- 1.0
P_min_tune <- 1.0
P_max_tune <- 1.0
gamma_T_tune <- 0.5
xi2_T_tune <- 2.5
gamma_P_tune <- 0.5
xi2_P_tune <- 2.85
W_T_tune <- 0.48
W_P_tune <- 0.28

T_obs_idx <- 447:556 
P_obs_idx <- 447:556

n_mcmc <- 15000
n_thin <- 5
psi <- 0.5

params <- list(n_mcmc = n_mcmc, num_species = num_species, mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0, mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1, alpha_s2 = alpha_s2, beta_s2 = beta_s2, alpha_Temp_min = alpha_Temp_min, beta_Temp_min = beta_Temp_min, Temp_min_lower = Temp_min_lower, Temp_min_upper = Temp_min_upper, alpha_Temp_max = alpha_Temp_max, beta_Temp_max = beta_Temp_max, Temp_max_lower = Temp_max_lower, Temp_max_upper = Temp_max_upper, alpha_P_min = alpha_P_min, beta_P_min = beta_P_min, P_min_lower = P_min_lower, P_min_upper = P_min_upper, alpha_P_max = alpha_P_max, beta_P_max = beta_P_max, P_max_lower = P_max_lower, P_max_upper = P_max_upper, mu_mu_gamma_T = mu_mu_gamma_T, s2_mu_gamma_T = s2_mu_gamma_T, alpha_s2_gamma_T = alpha_s2_gamma_T, beta_s2_gamma_T = beta_s2_gamma_T, mu_mu_xi2_T = mu_mu_xi2_T, s2_mu_xi2_T = s2_mu_xi2_T, alpha_s2_xi2_T = alpha_s2_xi2_T, beta_s2_xi2_T = beta_s2_xi2_T, mu_mu_gamma_P = mu_mu_gamma_P, s2_mu_gamma_P = s2_mu_gamma_P, alpha_s2_gamma_P = alpha_s2_gamma_P, beta_s2_gamma_P = beta_s2_gamma_P,  mu_mu_xi2_P = mu_mu_xi2_P, s2_mu_xi2_P = s2_mu_xi2_P, alpha_s2_xi2_P = alpha_s2_xi2_P, beta_s2_xi2_P = beta_s2_xi2_P, Temp_min_tune = Temp_min_tune, Temp_max_tune = Temp_max_tune, P_min_tune = P_min_tune, P_max_tune = P_max_tune, gamma_T_tune = gamma_T_tune, xi2_T_tune = xi2_T_tune, gamma_P_tune = gamma_P_tune, xi2_P_tune = xi2_P_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, day_len = day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, psi = psi, n_thin = n_thin)

##
## Run MCMC
##

out <- vector('list', length = 3)
for(i in 1:3){
	## Fit full model
  out[[i]] <- makeMCMC(y, Temp, P, species, params, process)
}

save(out, file = '~/Mechanistic-Tree-Ring/data/fitMixedHierarchicalRandomHudson.RData')
