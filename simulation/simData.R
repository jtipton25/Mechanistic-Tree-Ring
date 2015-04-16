simData <- function(model = 'VS', s2_noise = 0.5, trend = FALSE){
	##
	## Libraries and functions
	##
	
	library(myFunctions)
	library(bayesTreeRing)
	
	##
	## Load MCAR Data <- use parameters from the data and try to re-create
	##
	
	load('~/Mechanistic-Tree-Ring/datafiles/climateMCARTrend.RData')
	
	phi_vec <- c(rep(median(out$phi_one), 12), rep(median(out$phi_two), 12))
	Q_1_inv <- myFunctions::makeQinv(median(out$rho_one), 12)
	L_1 <- chol(Q_1_inv)
	Q_2_inv <- myFunctions::makeQinv(median(out$rho_two), 12)
	L_2 <- chol(Q_2_inv)
	L_1_L_2 <- t(L_1) %*% L_2
	Sigma_inv <- 1 / median(out$sigma_squared_W) * 
							rbind(cbind(Q_1_inv, median(out$omega) * L_1_L_2), 
							cbind(median(out$omega) * t(L_1_L_2), Q_2_inv))
	Sigma <- solve(Sigma_inv)
	det_Sigma <- - determinant(Sigma_inv, log = TRUE)$modulus[1]	
	
	##
	## Simulate the data using the observed data as a baseline
	##
	
	load('~/Mechanistic-Tree-Ring/datafiles/hudsonValleyData.RData')
	y <- as.matrix(crn.dat[ - 1, 2:35])
	H <- !is.na(y)
	
	year_idx <- rep(1:12, 110)
	
	Temp <- as.vector(t(matrix(base::rowMeans(
				Temp.avg.dat[ - which(Temp.avg.dat$Year > 2004), 3:36]), 110, 12)))
	mu_T <- as.vector(by(Temp, year_idx, mean))
	s_T <- as.vector(by(Temp, year_idx, sd))
	
	P <- as.vector(t(log(matrix(base::rowMeans(
				Precip.dat[ - which(Precip.dat$Year > 2004), 3:36]), 110, 12))))
	mu_P <- as.vector(by(P, year_idx, mean))
	s_P <- as.vector(by(P, year_idx, sd))
	
	day <- as.matrix(read.table('~/Mechanistic-Tree-Ring/datafiles/daylength', sep = "", 
										header = TRUE))
	day[day == -9999] <- NA
	day_len <- apply(day, 2, mean, na.rm = TRUE)
	day_len <- day_len / max(day_len)
	
	##
	## Initialize parameters
	##
	t <- 556
	obs <- 110
	num_species <- 12
	p <- 34
	W_tilde <- matrix(0, 24, t)
	W <- matrix(0, 24, t)
	W_tilde_0 <- rnorm(24, 0, 1)
	
	if(trend){ ## trend in the observation period not in reconstruction interval
		W_tilde[, 1] <- myFunctions::mvrnormArma(1, phi_vec * W_tilde_0, Sigma)
		W[1:12, 1] <- W_tilde[1:12, 1] * s_T + mu_T
		W[13:24, 1] <- W_tilde[13:24, 1] * s_P + mu_P
		for(k in 2:(t - obs)){
			W_tilde[, k] <- myFunctions::mvrnormArma(1, phi_vec * W_tilde[, k - 1], Sigma)
			W[1:12, k] <- W_tilde[1:12, k] * s_T + mu_T
			W[13:24, k] <- W_tilde[13:24, k] * s_P + mu_P
		}
		for(k in (t - obs + 1):t){
			W_tilde[, k] <- myFunctions::mvrnormArma(1, phi_vec * W_tilde[, k - 1], Sigma) + 
											c(rep( (k - (t - obs) + 1) / obs, 12), rep(0, 12))
			W[1:12, k] <- W_tilde[1:12, k] * s_T + mu_T
			W[13:24, k] <- W_tilde[13:24, k] * s_P + mu_P
		}
	} else { ## no trend in the simulation
		W_tilde[, 1] <- myFunctions::mvrnormArma(1, phi_vec * W_tilde_0, Sigma)
		W[1:12, 1] <- W_tilde[1:12, 1] * s_T + mu_T
		W[13:24, 1] <- W_tilde[13:24, 1] * s_P + mu_P
		for(k in 2:t){
			W_tilde[, k] <- myFunctions::mvrnormArma(1, phi_vec * W_tilde[, k - 1], Sigma)
			W[1:12, k] <- W_tilde[1:12, k] * s_T + mu_T
			W[13:24, k] <- W_tilde[13:24, k] * s_P + mu_P
		}
	}
	
	if (model == 'VS'){
		
		##
		## Set up VS-lite model parameters and simulate chronology given data
		##
		
		alpha_Temp_min = 9
		beta_Temp_min = 5
		Temp_min_lower = 0
		Temp_min_upper = 9
		alpha_Temp_max = 3.5
		beta_Temp_max = 3.5
		Temp_max_lower = 10
		Temp_max_upper = 24
		alpha_P_min = 3.5
		beta_P_min = 3.5
		P_min_lower = 65
		P_min_upper = 85
		alpha_P_max = 3.5
		beta_P_max = 3.5
		P_max_lower = 85
		P_max_upper = 105
		
		species <- c(1:num_species, sample(1:num_species, p-num_species, replace = TRUE))
		species_minus_one <- species - 1
		
		Temp_min <- rbeta(num_species, alpha_Temp_min, beta_Temp_min) * 
								(Temp_min_upper - Temp_min_lower) + Temp_min_lower
		Temp_max <- rbeta(num_species, alpha_Temp_max, beta_Temp_max) * 
								(Temp_max_upper - Temp_max_lower) + Temp_max_lower
		P_min <- rbeta(num_species, alpha_P_min, beta_P_max) * (P_min_upper - P_min_lower) + 
							P_min_lower
		P_max <- rbeta(num_species, alpha_P_max, beta_P_max) * (P_max_upper - P_max_lower) + 
							P_max_lower
		
		zeta_pre <- bayesTreeRing::makeZetaVS(t, p, day_len, W, Temp_min, Temp_max, P_min, P_max, 
												 species_minus_one)
		zeta <- matrix(0, t, p)
		for(k in 1:p){
			zeta_tmp <- zeta_pre[, k]
			zeta_mn <- mean(zeta_tmp)
			zeta_sd <- sd(zeta_tmp)
			zeta[, k]	<- 1 / 5 * rnorm(t, sqrt(1 - s2_noise) * (zeta_tmp - zeta_mn) / zeta_sd, sqrt(s2_noise)) + 1
		}
		list(zeta = zeta, Temp_min = Temp_min, Temp_max = Temp_max, P_min = P_min, 
				 P_max = P_max, H = H, t = t, p = p, num_species = num_species, W = W, W_tilde = W_tilde, 
				 species = species, day_len = day_len)
		
	} else if (model == 'Probit'){
		
		##
		## Set up Probit model parameters and simulate chronology given data
		##
		
		mu_gamma_T = 12
		sigma_squared_gamma_T = 8
		mu_xi_squared_T = log(4)
		sigma_squared_xi_squared_T = 0.125
		mu_gamma_P = 85
		sigma_squared_gamma_P = 64
		mu_xi_squared_P = log(6)
		sigma_squared_xi_squared_P = 0.125
		
		species <- c(1:num_species, sample(1:num_species, p-num_species, replace = TRUE))
		species_minus_one <- species - 1
		
		gamma_T <- rnorm(num_species, mu_gamma_T, sqrt(sigma_squared_gamma_T))
		xi_squared_T <- rlnorm(num_species, mu_xi_squared_T, sqrt(sigma_squared_xi_squared_T))
		xi_T <- sqrt(xi_squared_T)
		gamma_P <- rnorm(num_species, mu_gamma_P, sqrt(sigma_squared_gamma_P))
		xi_squared_P <- rlnorm(num_species, mu_xi_squared_P, sqrt(sigma_squared_xi_squared_P))
		xi_P <- sqrt(xi_squared_P)
		
		zeta_pre <- bayesTreeRing::makeZetaPro(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
												 species_minus_one)
		
		zeta <- matrix(0, t, p)
		for(k in 1:p){
			zeta_tmp <- zeta_pre[, k]
			zeta_mn <- mean(zeta_tmp)
			zeta_sd <- sd(zeta_tmp)
			zeta[, k]	<- 1 / 5 * rnorm(t, sqrt(1 - s2_noise) * (zeta_tmp - zeta_mn) / zeta_sd, sqrt(s2_noise)) + 1
		}
		
		list(zeta = zeta, gamma_T = gamma_T, xi_squared_T = xi_squared_T, 
				 gamma_P = gamma_P, xi_squared_P = xi_squared_P, H = H, t = t, p = p, 
				 num_species = num_species, W = W, W_tilde = W_tilde, species = species, day_len = day_len)
		
		
	} else if (model == 'Mixed'){
	##
	## Simulate mixed VS-Lite and Probit chronology data
	##
	
	## Set up VS-Lite model parameters 
	alpha_Temp_min = 9
	beta_Temp_min = 5
	Temp_min_lower = 0
	Temp_min_upper = 9
	alpha_Temp_max = 3.5
	beta_Temp_max = 3.5
	Temp_max_lower = 10
	Temp_max_upper = 24
	alpha_P_min = 3.5
	beta_P_min = 3.5
	P_min_lower = 65
	P_min_upper = 85
	alpha_P_max = 3.5
	beta_P_max = 3.5
	P_max_lower = 85
	P_max_upper = 105
	## Set up Probit model parameters 
	mu_gamma_T = 12
	sigma_squared_gamma_T = 8
	mu_xi_squared_T = log(4)
	sigma_squared_xi_squared_T = 0.125
	mu_gamma_P = 85
	sigma_squared_gamma_P = 64
	mu_xi_squared_P = log(6)
	sigma_squared_xi_squared_P = 0.125
	
	species <- c(1:num_species, sample(1:num_species, p-num_species, replace = TRUE))
	species_minus_one <- species - 1
	
	Temp_min <- rbeta(num_species, alpha_Temp_min, beta_Temp_min) * 
							(Temp_min_upper - Temp_min_lower) + Temp_min_lower
	Temp_max <- rbeta(num_species, alpha_Temp_max, beta_Temp_max) * 
							(Temp_max_upper - Temp_max_lower) + Temp_max_lower
	P_min <- rbeta(num_species, alpha_P_min, beta_P_max) * (P_min_upper - P_min_lower) + 
								P_min_lower
	P_max <- rbeta(num_species, alpha_P_max, beta_P_max) * (P_max_upper - P_max_lower) +
								P_max_lower
	
	gamma_T <- rnorm(num_species, mu_gamma_T, sqrt(sigma_squared_gamma_T))
	xi_squared_T <- rlnorm(num_species, mu_xi_squared_T, sqrt(sigma_squared_xi_squared_T))
	xi_T <- sqrt(xi_squared_T)
	gamma_P <- rnorm(num_species, mu_gamma_P, sqrt(sigma_squared_gamma_P))
	xi_squared_P <- rlnorm(num_species, mu_xi_squared_P, sqrt(sigma_squared_xi_squared_P))
	xi_P <- sqrt(xi_squared_P)
	x <- matrix(rbinom(t * num_species, 1, 0.5), t, num_species)
	zeta_pre <- bayesTreeRing::makeZetaMix(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, Temp_min, 
											 Temp_max, P_min, P_max, species_minus_one, x)
	zeta <- matrix(0, t, p)
	for(k in 1:p){
		zeta_tmp <- zeta_pre[, k]
		zeta_mn <- mean(zeta_tmp)
		zeta_sd <- sd(zeta_tmp)
		zeta[, k]	<- 1 / 5 * rnorm(t, sqrt(1 - s2_noise) * (zeta_tmp - zeta_mn) / zeta_sd, sqrt(s2_noise)) + 1
	}
	list(zeta = zeta, gamma_T = gamma_T, xi_squared_T = xi_squared_T,
			 gamma_P = gamma_P, xi_squared_P = xi_squared_P, Temp_min = Temp_min,
			 Temp_max = Temp_max, P_min = P_min, P_max = P_max, H = H, t = t, p = p, 
       num_species = num_species, W = W, W_tilde = W_tilde, species = species, day_len = day_len, 
			 x = x)
	} else if(model == 'normalProbit'){
		
	  
	  ##
	  ## Set up model parameters and simulate chronology given data
	  ##
	  
	  mu_gamma_T = 16
	  sigma_squared_gamma_T = 8
	  mu_xi_squared_T = log(4)
	  sigma_squared_xi_squared_T = 0.125
	  mu_gamma_P = 85
	  sigma_squared_gamma_P = 64
	  mu_xi_squared_P = log(6)
	  sigma_squared_xi_squared_P = 0.125
	  
	  species <- c(1:num_species, sample(1:num_species, p-num_species, replace = TRUE))
	  species_minus_one <- species - 1
	  
	  gamma_T <- rnorm(num_species, mu_gamma_T, sqrt(sigma_squared_gamma_T))
	  xi_squared_T <- rlnorm(num_species, mu_xi_squared_T, sqrt(sigma_squared_xi_squared_T))
	  xi_T <- sqrt(xi_squared_T)
	  gamma_P <- rnorm(num_species, mu_gamma_P, sqrt(sigma_squared_gamma_P))
	  xi_squared_P <- rlnorm(num_species, mu_xi_squared_P, sqrt(sigma_squared_xi_squared_P))
	  xi_P <- sqrt(xi_squared_P)
	  
	  zeta_pre <- bayesTreeRing::makeZetaNP(t, p, day_len, W, gamma_T, xi_T, gamma_P, xi_P, 
	                       species_minus_one)
	  
	  zeta <- matrix(0, t, p)
	  for(k in 1:p){
	    zeta_tmp <- zeta_pre[, k]
	    zeta_mn <- mean(zeta_tmp)
	    zeta_sd <- sd(zeta_tmp)
	    zeta[, k]	<- 1 / 5 * rnorm(t, sqrt(1 - s2_noise) * (zeta_tmp - zeta_mn) / zeta_sd, sqrt(s2_noise)) + 1
	  }
	  
	  list(zeta = zeta, gamma_T = gamma_T, xi_squared_T = xi_squared_T, 
	       gamma_P = gamma_P, xi_squared_P = xi_squared_P, H = H, t = t, p = p, 
	       num_species = num_species, W = W, W_tilde = W_tilde, species = species, day_len = day_len)
	}
}
