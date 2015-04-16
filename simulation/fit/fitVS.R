makeVSHierarchical <- function(dataSim, model, process, trend, n_mcmc){
	
	##
	## Libraries and functions
	##
	
	Rcpp::sourceCpp('~/Mechanistic-Tree-Ring/mcmc/mcmcVS.cpp')
	
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
	
	##VS_Lite growth model parameters
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
	
	## tuning parameters
	Temp_min_tune <- 0.75
	Temp_max_tune <- 0.75
	P_min_tune <- 1.4
	P_max_tune <- 1.75
	W_T_tune <- 0.46
	W_P_tune <- 0.25
	
	T_obs_idx <- 447:556 
	P_obs_idx <- 447:556
	
	n_thin <- 5
	
	params <- list(n_mcmc = n_mcmc, num_species = dataSim$num_species, mu_mu_beta_0 = mu_mu_beta_0, s2_mu_beta_0 = s2_mu_beta_0, alpha_s2_beta_0 = alpha_s2_beta_0, beta_s2_beta_0 = beta_s2_beta_0, mu_mu_beta_1 = mu_mu_beta_1, s2_mu_beta_1 = s2_mu_beta_1, alpha_s2_beta_1 = alpha_s2_beta_1, beta_s2_beta_1 = beta_s2_beta_1, alpha_s2 = alpha_s2, beta_s2 = beta_s2, alpha_Temp_min = alpha_Temp_min, beta_Temp_min = beta_Temp_min, alpha_Temp_max = alpha_Temp_max, beta_Temp_max = beta_Temp_max, alpha_P_min = alpha_P_min, beta_P_min = beta_P_min, alpha_P_max = alpha_P_max, beta_P_max = beta_P_max, Temp_min_lower = Temp_min_lower, Temp_min_upper = Temp_min_upper, Temp_max_lower = Temp_max_lower, Temp_max_upper = Temp_max_upper, P_min_lower = P_min_lower, P_min_upper = P_min_upper, P_max_lower = P_max_lower, P_max_upper = P_max_upper, Temp_min_tune = Temp_min_tune, Temp_max_tune = Temp_max_tune, P_min_tune = P_min_tune,  P_max_tune = P_max_tune, W_T_tune = W_T_tune, W_P_tune = W_P_tune, day_len = dataSim$day_len, T_obs_idx = T_obs_idx, P_obs_idx = P_obs_idx, n_thin = n_thin)
	
	
	##
	## Run MCMC
	##
	
	## Fit models
	
	start <- Sys.time()
	if(model == 'VS'){
		simVSfitVSHierarchical <- vector('list', length = 3)
		for(k  in 1:3){
			simVSfitVSHierarchical[[k]] <- makeMCMC(dataSim$zeta * dataSim$H,
																							dataSim$W[1:12, 447:556],
																							dataSim$W[13:24, 447:556], dataSim$species,
																							params, process, sim = TRUE)
		}
		save(simVSfitVSHierarchical, dataSim, file = 
				 	'~/Mechanistic-Tree-Ring/data/simVSfitVSHierarchical.RData')
	} else if (model == 'Probit'){
		simProbitfitVSHierarchical <- vector('list', length = 3)
		for(k  in 1:3){
			simProbitfitVSHierarchical[[k]] <- makeMCMC(dataSim$zeta * dataSim$H,
																									dataSim$W[1:12, 447:556], 
																									dataSim$W[13:24, 447:556], dataSim$species, 
																									params, process, sim = TRUE)
		}
		save(simProbitfitVSHierarchical, dataSim, file = 
				 	'~/Mechanistic-Tree-Ring/simProbitfitVSHierarchical.RData')
	} else if (model == 'normalProbit'){
		simNormalProbitfitVSHierarchical <- vector('list', length = 3)
		for(k in 1:3){
			simNormalProbitfitVSHierarchical[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556], dataSim$species, params, process, sim = TRUE)
		}
		save(simNormalProbitfitVSHierarchical, dataSim, file = '~/Mechanistic-Tree-Ring/data/simNormalProbitfitVSHierarchical.RData')
	} else {
		simMixedfitVSHierarchical <- vector('list', length = 3)
		for(k  in 1:3){
			simMixedfitVSHierarchical[[k]] <- makeMCMC(dataSim$zeta * dataSim$H, 
																								 dataSim$W[1:12, 447:556],
																								 dataSim$W[13:24, 447:556], dataSim$species, 
																								 params, process, sim = TRUE)
		}
		save(simMixedfitVSHierarchical, dataSim, 
				 file = '~/Mechanistic-Tree-Ring/data/simMixedfitVSHierarchical.RData')
	}
	finish <- Sys.time()
	finish - start
}

