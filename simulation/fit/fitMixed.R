makeMixedHierarchicalRandom <- function(dataSim, model, process, trend, n_mcmc){
	
	##
	## Libraries and functions
	##
	
	Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcMixed.cpp')
	
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
	
	n_thin <- 5
	psi <- 0.5
	
	params <- list(n_mcmc = n_mcmc, num_species = dataSim$num_species, mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0, mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1, alpha_s2 = alpha_s2, beta_s2 = beta_s2, alpha_Temp_min = alpha_Temp_min, beta_Temp_min = beta_Temp_min, Temp_min_lower = Temp_min_lower, Temp_min_upper = Temp_min_upper, alpha_Temp_max = alpha_Temp_max, beta_Temp_max = beta_Temp_max, Temp_max_lower = Temp_max_lower, Temp_max_upper = Temp_max_upper, alpha_P_min = alpha_P_min, beta_P_min = beta_P_min, P_min_lower = P_min_lower, P_min_upper = P_min_upper, alpha_P_max = alpha_P_max, beta_P_max = beta_P_max, P_max_lower = P_max_lower, P_max_upper = P_max_upper, mu_mu_gamma_T = mu_mu_gamma_T, s2_mu_gamma_T = s2_mu_gamma_T, alpha_s2_gamma_T = alpha_s2_gamma_T, beta_s2_gamma_T = beta_s2_gamma_T, mu_mu_xi2_T = mu_mu_xi2_T, s2_mu_xi2_T = s2_mu_xi2_T, alpha_s2_xi2_T = alpha_s2_xi2_T, beta_s2_xi2_T = beta_s2_xi2_T, mu_mu_gamma_P = mu_mu_gamma_P, s2_mu_gamma_P = s2_mu_gamma_P, alpha_s2_gamma_P = alpha_s2_gamma_P, beta_s2_gamma_P = beta_s2_gamma_P,  mu_mu_xi2_P = mu_mu_xi2_P, s2_mu_xi2_P = s2_mu_xi2_P, alpha_s2_xi2_P = alpha_s2_xi2_P, beta_s2_xi2_P = beta_s2_xi2_P, Temp_min_tune = Temp_min_tune, Temp_max_tune = Temp_max_tune, P_min_tune = P_min_tune, P_max_tune = P_max_tune, gamma_T_tune = gamma_T_tune, xi2_T_tune = xi2_T_tune, gamma_P_tune = gamma_P_tune, xi2_P_tune = xi2_P_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, day_len = dataSim$day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, psi = psi, n_thin = n_thin)
	
	
	##
	## Run MCMC
	##
	
	## model fit
	start <- Sys.time()
	if (model == 'VS'){
		simVSfitMixedHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simVSfitMixedHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simVSfitMixedHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simVSfitMixedHierarchicalRandom.RData')
	} else if (model == 'Probit'){
		simProbitfitMixedHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simProbitfitMixedHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simProbitfitMixedHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simProbitfitMixedHierarchicalRandom.RData')
	} else if (model == 'normalProbit'){
		simNormalProbitfitMixedHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simNormalProbitfitMixedHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simNormalProbitfitMixedHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simNormalProbitfitMixedHierarchicalRandom.RData')
	} else {
		simMixedfitMixedHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simMixedfitMixedHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simMixedfitMixedHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simMixedfitMixedHierarchicalRandom.RData')
	}
	
	finish <- Sys.time()
	finish - start
	
	rm(list = ls())
}
