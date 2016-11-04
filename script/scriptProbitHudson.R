set.seed(111)
##
## Libraries and functions
##

Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcProbit.cpp')

library(bayesTreeRing)
library(myFunctions)

##
## load the data
##
data("hudsonValleyData")
# load('~/Mechanistic-Tree-Ring/datafiles/hudsonValleyData.RData')
species <- as.integer(as.factor(species.crn[ - 1]))
num_species <- length(unique(species))
p <- length(species) 
y <- as.matrix(crn.dat[ - 1, 2:35])
H <- !is.na(y)
y[!H] <- 0L


J <- rep(1, 110)
year_idx <- rep(1:12, 110)

Temp <- t(matrix(base::rowMeans(
	Temp.avg.dat[ - which(Temp.avg.dat$Year > 2004), 3:36]), 110, 12))
P <- t(log(matrix(base::rowMeans(
	Precip.dat[ - which(Precip.dat$Year > 2004), 3:36]), 110, 12)))

day <- as.matrix(read.table('~/Mechanistic-Tree-Ring/datafiles/daylength', 
														sep = "", header = TRUE))
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
mu_mu_gamma_T <- 13
s2_mu_gamma_T <- 100
# s2_mu_gamma_T <- 16
alpha_s2_gamma_T <- 1
beta_s2_gamma_T <- 1
# alpha_s2_gamma_T <- 5
# beta_s2_gamma_T <- 0.2
mu_mu_xi2_T <- 0
s2_mu_xi2_T <- 10
# s2_mu_xi2_T <- 0.5
alpha_s2_xi2_T <- 1
beta_s2_xi2_T <- 0.1
# alpha_s2_xi2_T <- 10
# beta_s2_xi2_T <- 0.1
mu_mu_gamma_P <- 85
s2_mu_gamma_P <- 36
alpha_s2_gamma_P <- 1
beta_s2_gamma_P <- 1
# alpha_s2_gamma_P <- 5
# beta_s2_gamma_P <- 0.2
mu_mu_xi2_P <- 0
s2_mu_xi2_P <- 1
alpha_s2_xi2_P <- 1
beta_s2_xi2_P <- 1
# alpha_s2_xi2_P <- 10
# beta_s2_xi2_P <- 0.1

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
  
params <- list(n_mcmc = n_mcmc, num_species = num_species, 
							 mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, 
							 alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0,
							 mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, 
							 alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1,
							 alpha_s2 = alpha_s2, beta_s2 = beta_s2, mu_mu_gamma_T = mu_mu_gamma_T,
							 s2_mu_gamma_T = s2_mu_gamma_T, alpha_s2_gamma_T = alpha_s2_gamma_T, 
							 beta_s2_gamma_T = beta_s2_gamma_T, mu_mu_xi2_T = mu_mu_xi2_T, 
							 s2_mu_xi2_T = s2_mu_xi2_T, alpha_s2_xi2_T = alpha_s2_xi2_T,
							 beta_s2_xi2_T = beta_s2_xi2_T, mu_mu_gamma_P = mu_mu_gamma_P, 
							 s2_mu_gamma_P = s2_mu_gamma_P, alpha_s2_gamma_P = alpha_s2_gamma_P,
							 beta_s2_gamma_P = beta_s2_gamma_P, mu_mu_xi2_P = mu_mu_xi2_P, 
							 s2_mu_xi2_P = s2_mu_xi2_P, alpha_s2_xi2_P = alpha_s2_xi2_P, 
							 beta_s2_xi2_P = beta_s2_xi2_P, gamma_T_tune = gamma_T_tune, 
							 xi2_T_tune = xi2_T_tune, gamma_P_tune = gamma_P_tune, 
							 xi2_P_tune = xi2_P_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, 
							 day_len = day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, 
							 n_thin = n_thin)

##
## Run MCMC
##

start <- Sys.time()
out <- makeMCMC(y, Temp, P, species, params, process)
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
matplot(apply(W_T_month, 1, median), type = 'l', ylim = c(6, 12),
				main = round(out$accept$W_T_accept, 4), ylab = 'T', xaxt = 'n',
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_T_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_T_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)

W_P_month <- downscaleMonth(out$W[13:24, , ])
matplot(apply(W_P_month, 1, median), type = 'l', ylim = c(4, 5),
				main = round(out$accept$W_P_accept, 4), ylab = 'log(P)', xaxt = 'n', 
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_P_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_P_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
pederson <- read.table(file = '~/Google Drive/NYCDroughtPederson.txt', 
											 sep = '', header = TRUE)
tmp <- rbind(matrix(NA, 82, 2, dimnames = list(c(), c("YEAR", "drought"))),
						 pederson[1:474, ])
matplot((tmp$drought / sd(pederson$drought)) * sd(W_P_month) + mean(W_P_month),
				type = 'l', col = 'blue', add = TRUE)

##
## Growing Season plots
##

layout(matrix(1:2, 2))

W_T_summer_month <- downscaleMonth(out$W[4:10, , ])
matplot(apply(W_T_summer_month, 1, median), type = 'l', ylim = c(14, 18),
				main = round(out$accept$W_T_accept, 4), ylab = 'T', xaxt = 'n',
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_T_summer_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_T_summer_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)

W_P_summer_month <- downscaleMonth(out$W[16:22, , ])
matplot(apply(W_P_summer_month, 1, median), type = 'l', ylim = c(4, 5),
				main = round(out$accept$W_P_accept, 4), ylab = 'log(P)', xaxt = 'n', 
				xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(W_P_summer_month, 1, quantile, prob = 0.975), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(W_P_summer_month, 1, quantile, prob = 0.025), type = 'l', 
				col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot((tmp$drought / sd(pederson$drought)) * sd(W_P_summer_month) + 
					mean(W_P_summer_month), type = 'l', col = 'blue', add = TRUE)


