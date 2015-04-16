makeProbitHierarchicalRandom <- function(dataSim, model, process, trend, n_mcmc){
	
	##
	## Libraries and functions
	##
	
	Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcProbit.cpp')
	
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
	
	n_thin <- 5
	
	params <- list(n_mcmc = n_mcmc, num_species = dataSim$num_species, mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0, mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1, alpha_s2 = alpha_s2, beta_s2 = beta_s2, mu_mu_gamma_T = mu_mu_gamma_T, s2_mu_gamma_T = s2_mu_gamma_T, alpha_s2_gamma_T = alpha_s2_gamma_T, beta_s2_gamma_T = beta_s2_gamma_T, mu_mu_xi2_T = mu_mu_xi2_T, s2_mu_xi2_T = s2_mu_xi2_T, alpha_s2_xi2_T = alpha_s2_xi2_T, beta_s2_xi2_T = beta_s2_xi2_T, mu_mu_gamma_P = mu_mu_gamma_P, s2_mu_gamma_P = s2_mu_gamma_P, alpha_s2_gamma_P = alpha_s2_gamma_P, beta_s2_gamma_P = beta_s2_gamma_P, mu_mu_xi2_P = mu_mu_xi2_P, s2_mu_xi2_P = s2_mu_xi2_P, alpha_s2_xi2_P = alpha_s2_xi2_P, beta_s2_xi2_P = beta_s2_xi2_P, gamma_T_tune = gamma_T_tune, xi2_T_tune = xi2_T_tune, gamma_P_tune = gamma_P_tune, xi2_P_tune = xi2_P_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, day_len = dataSim$day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, n_thin = n_thin)
	
	
	
	##
	## Run MCMC
	##
	## Fit model
	
	start <- Sys.time()
	if(model == 'VS'){
		simVSfitProbitHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simVSfitProbitHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, 
																			dataSim$W[1:12, 447:556], 
																			dataSim$W[13:24, 447:556], 
																			dataSim$species, params, process, 
																			sim = TRUE)
		}
		save(simVSfitProbitHierarchicalRandom, dataSim, file = 
				 	'~/Mechanistic-Tree-Ring/data/simVSfitProbitHierarchicalRandom.RData')
	} else if (model == 'Probit'){
		simProbitfitProbitHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simProbitfitProbitHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, 
																					dataSim$W[1:12, 447:556], 
																					dataSim$W[13:24, 447:556], 
																					dataSim$species, params, process, 
																					sim = TRUE)
		}
		save(simProbitfitProbitHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simProbitfitProbitHierarchicalRandom.RData')
	} else if (model == 'normalProbit'){
		simNormalProbitfitProbitHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simNormalProbitfitProbitHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simNormalProbitfitProbitHierarchicalRandom, dataSim, file = '~/Mechanistic-Tree-Ring/data/simNormalProbitfitProbitHierarchicalRandom.RData')
	} else {
		simMixedfitProbitHierarchicalRandom <- vector('list', length = 3)
		for(k in 1:3){
			simMixedfitProbitHierarchicalRandom[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, 
																				 dataSim$W[1:12, 447:556],
																				 dataSim$W[13:24, 447:556], 
																				 dataSim$species, params, process, 
																				 sim = TRUE)
		}
		save(simMixedfitProbitHierarchicalRandom, dataSim, file =
				 	'~/Mechanistic-Tree-Ring/data/simMixedfitProbitHierarchicalRandom.RData')
	}
	finish <- Sys.time()
	finish - start
	rm(list = ls())
}
