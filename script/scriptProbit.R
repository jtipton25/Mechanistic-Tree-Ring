set.seed(101)
##
## Libraries and functions
##

Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcProbit.cpp')

library(bayesTreeRing)
library(myFunctions)

## simulate data
dataSim <- simData(model = 'Probit', s2_noise = 0.75, trend = TRUE)

## fit empirical Bayes climate model
process <- makeMCARTrend(dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556])

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
gamma_T_tune <- 0.5
xi2_T_tune <- 3.5
gamma_P_tune <- 0.75
xi2_P_tune <- 3.0
W_T_tune <- 0.40
W_P_tune <- 0.18

T_obs_idx <- 447:556 
P_obs_idx <- 447:556
  
n_mcmc <- 5000
n_thin <- 5
  
params <- list(n_mcmc = n_mcmc, num_species = dataSim$num_species, mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0, mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1, alpha_s2 = alpha_s2, beta_s2 = beta_s2, mu_mu_gamma_T = mu_mu_gamma_T, s2_mu_gamma_T = s2_mu_gamma_T, alpha_s2_gamma_T = alpha_s2_gamma_T, beta_s2_gamma_T = beta_s2_gamma_T, mu_mu_xi2_T = mu_mu_xi2_T, s2_mu_xi2_T = s2_mu_xi2_T, alpha_s2_xi2_T = alpha_s2_xi2_T, beta_s2_xi2_T = beta_s2_xi2_T, mu_mu_gamma_P = mu_mu_gamma_P, s2_mu_gamma_P = s2_mu_gamma_P, alpha_s2_gamma_P = alpha_s2_gamma_P, beta_s2_gamma_P = beta_s2_gamma_P, mu_mu_xi2_P = mu_mu_xi2_P, s2_mu_xi2_P = s2_mu_xi2_P, alpha_s2_xi2_P = alpha_s2_xi2_P, beta_s2_xi2_P = beta_s2_xi2_P, gamma_T_tune = gamma_T_tune, xi2_T_tune = xi2_T_tune, gamma_P_tune = gamma_P_tune, xi2_P_tune = xi2_P_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, day_len = dataSim$day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, n_thin = n_thin)

##
## Run MCMC
##

start <- Sys.time()
out <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, 
								params, process, sim = TRUE)
finish <- Sys.time()
finish - start

##
## Plot MCMC
##

layout(matrix(1:6, 3, 2))
matplot(t(out$calibration$beta_0), type = 'l')
matplot(out$calibration$mu_beta_0, type = 'l')
matplot(out$calibration$s2_beta_0, type = 'l')
matplot(t(out$calibration$beta_1), type = 'l')
matplot(out$calibration$mu_beta_1, type = 'l')
matplot(out$calibration$s2_beta_1, type = 'l')

layout(matrix(1))
matplot(t(out$calibration$s2), type = 'l')

layout(matrix(1:6, 3, 2)) 
matplot(t(out$gamma_T), type = 'l', 
				main = round(mean(out$accept$gamma_T_accept), 4))
matplot(out$mu_gamma_T, type = 'l')
matplot(out$s2_gamma_T, type = 'l')
matplot(t(sqrt(out$xi2_T)), type = 'l', 
				main = round(mean(out$accept$xi2_T_accept), 4))
matplot(out$mu_xi2_T, type = 'l')
matplot(out$s2_xi2_T, type = 'l')

layout(matrix(1:6, 3, 2)) 
matplot(t(out$gamma_P), type = 'l', 
				main = round(mean(out$accept$gamma_P_accept), 4))
matplot(out$mu_gamma_P, type = 'l')
matplot(out$s2_gamma_P, type = 'l')
matplot(t(sqrt(out$xi2_P)), type = 'l', 
				main = round(mean(out$accept$xi2_P_accept), 4))
matplot(out$mu_xi2_P, type = 'l')
matplot(out$s2_xi2_P, type = 'l')

year <- 2005 - (556:1)
year_idx <- rep(1:556, each = 12)

layout(matrix(1:2, 2))
W_T_month <- downscaleMonth(out$W[1:12, , ])
W_T_true <- apply(dataSim$W[1:12, ], 2, mean)
matplot(apply(W_T_month, 1, median), type = 'l', ylim = c(6, 12),
				main = round(out$accept$W_T_accept, 4), ylab = 'T', xaxt = 'n', xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_T_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_T_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(W_T_true, type = 'l', add = TRUE, col = 'blue')

W_P_month <- downscaleMonth(out$W[13:24, , ])
W_P_true <- apply(dataSim$W[13:24, ], 2, mean)
matplot(apply(W_P_month, 1, median), type = 'l', ylim = c(4, 5),
				main = round(out$accept$W_P_accept, 4), ylab = 'log(P)', xaxt = 'n', 
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_P_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_P_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(W_P_true, type = 'l', add = TRUE, col = 'blue')

##
## Growing Season plots
##

layout(matrix(1:2, 2))

W_T_summer_month <- downscaleMonth(out$W[4:10, , ])
W_T_summer_true <- apply(dataSim$W[4:10, ], 2, mean)
matplot(apply(W_T_summer_month, 1, median), type = 'l', ylim = c(14, 18),
				main = round(out$accept$W_T_accept, 4), ylab = 'T', xaxt = 'n', xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_T_summer_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_T_summer_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(W_T_summer_true, type = 'l', add = TRUE, col = 'blue')

W_P_summer_month <- downscaleMonth(out$W[16:22, , ])
W_P_summer_true <- apply(dataSim$W[16:22, ], 2, mean)
matplot(apply(W_P_summer_month, 1, median), type = 'l', ylim = c(4, 5),
				main = round(out$accept$W_P_accept, 4), ylab = 'log(P)', xaxt = 'n', 
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_P_summer_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_P_summer_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(W_P_summer_true, type = 'l', add = TRUE, col = 'blue')

##
## CRPS score
##

## Climatology predictions
Temp_hat <- predictClimatology(n_mcmc / 2, dataSim$W[1:12, 447:556])
P_hat <- predictClimatology(n_mcmc / 2, dataSim$W[13:24, 447:556])

## CRPS for T relative to climatology
climate_CRPS_T <- makeCRPS(445, Temp_hat[ - 1, ], W_T_true[ -1 ])
mean(climate_CRPS_T)
model_CRPS_T <- makeCRPS(445, W_T_month[ - 1, ], W_T_true[ -1 ])
mean(model_CRPS_T)

## CRPS for P relative to climatology
climate_CRPS_P <- makeCRPS(445, P_hat[ - 1, ], W_P_true[ -1 ])
model_CRPS_P <- makeCRPS(445, W_P_month[ - 1, ], W_P_true[ -1 ])

library(ggplot2)
library(gridExtra)
CRPS_T = data.frame(CRPS = c(climate_CRPS_T, model_CRPS_T),
										Predictor = factor(rep(c("climatology", "model"), each = 445)))
CRPS_P = data.frame(CRPS = c(climate_CRPS_P, model_CRPS_P), 
										Predictor = factor(rep(c("climatology", "model"), each = 445)))

ggroup_T = ggplot(data = CRPS_T, aes(x = CRPS, fill = Predictor)) +
	geom_density(alpha = 0.3) + ggtitle("CRPS for T") + 
	guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
	theme(legend.key = element_rect(fill = "black"))
ggroup_P = ggplot(data = CRPS_P, aes(x = CRPS, fill = Predictor)) +
	geom_density(alpha = 0.3) + ggtitle("CRPS for P") + 
	guides(fill = guide_legend(override.aes = list(colour = NULL)))

multiplot(ggroup_T, ggroup_P, cols = 2)

