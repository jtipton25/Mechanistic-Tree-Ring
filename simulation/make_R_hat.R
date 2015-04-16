make_R_hat <- function(out, model, hierarchical = FALSE, random = FALSE){
	
	## R_hat_beta_0
	m <- length(out)
	n <- dim(out[[1]]$calibration$beta_0)[2]
	Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$calibration$beta_0), rowMeans(out[[2]]$calibration$beta_0), rowMeans(out[[3]]$calibration$beta_0)))
	B <- n / (m - 1) *((rowMeans(out[[1]]$calibration$beta_0) - Psi_bar)^2 + (rowMeans(out[[2]]$calibration$beta_0) - Psi_bar)^2 + (rowMeans(out[[3]]$calibration$beta_0) - Psi_bar)^2)
	W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$calibration$beta_0 - rowMeans(out[[1]]$calibration$beta_0))^2) + 1 / (n - 1) * rowSums((out[[2]]$calibration$beta_0 - rowMeans(out[[2]]$calibration$beta_0))^2) + 1 / (n - 1) * rowSums((out[[3]]$calibration$beta_0 - rowMeans(out[[3]]$calibration$beta_0))^2))
	R_hat_beta_0 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
	
	## R_hat_beta_1
	Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$calibration$beta_1), rowMeans(out[[2]]$calibration$beta_1), rowMeans(out[[3]]$calibration$beta_1)))
	B <- n / (m - 1) *((rowMeans(out[[1]]$calibration$beta_1) - Psi_bar)^2 + (rowMeans(out[[2]]$calibration$beta_1) - Psi_bar)^2 + (rowMeans(out[[3]]$calibration$beta_1) - Psi_bar)^2)
	W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$calibration$beta_1 - rowMeans(out[[1]]$calibration$beta_1))^2) + 1 / (n - 1) * rowSums((out[[2]]$calibration$beta_1 - rowMeans(out[[2]]$calibration$beta_1))^2) + 1 / (n - 1) * rowSums((out[[3]]$calibration$beta_1 - rowMeans(out[[3]]$calibration$beta_1))^2))
	R_hat_beta_1 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
	
	## R_hat for s2
	Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$calibration$s2), rowMeans(out[[2]]$calibration$s2), rowMeans(out[[3]]$calibration$s2)))
	B <- n / (m - 1) * ((rowMeans(out[[1]]$calibration$s2) - Psi_bar)^2 + (rowMeans(out[[2]]$calibration$s2) - Psi_bar)^2 + (rowMeans(out[[3]]$calibration$s2) - Psi_bar)^2)
	W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$calibration$s2 - rowMeans(out[[1]]$calibration$s2))^2) + 1 / (n - 1) * rowSums((out[[2]]$calibration$s2 - rowMeans(out[[2]]$calibration$s2))^2) + 1 / (n - 1) * rowSums((out[[3]]$calibration$s2 - rowMeans(out[[3]]$calibration$s2))^2))
	R_hat_s2 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
	
	if(hierarchical){
		## R_hat_mu_beta_0
		Psi_bar <- mean(cbind(mean(out[[1]]$calibration$mu_beta_0), mean(out[[2]]$calibration$mu_beta_0), mean(out[[3]]$calibration$mu_beta_0)))
		B <- n / (m - 1) *((mean(out[[1]]$calibration$mu_beta_0) - Psi_bar)^2 + (mean(out[[2]]$calibration$mu_beta_0) - Psi_bar)^2 + (mean(out[[3]]$calibration$mu_beta_0) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$calibration$mu_beta_0 - mean(out[[1]]$calibration$mu_beta_0))^2) + 1 / (n - 1) * sum((out[[2]]$calibration$mu_beta_0 - mean(out[[2]]$calibration$mu_beta_0))^2) + 1 / (n - 1) * sum((out[[3]]$calibration$mu_beta_0 - mean(out[[3]]$calibration$mu_beta_0))^2))
		R_hat_beta_0 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_s2_beta_0
		Psi_bar <- mean(cbind(mean(out[[1]]$calibration$s2_beta_0), mean(out[[2]]$calibration$s2_beta_0), mean(out[[3]]$calibration$s2_beta_0)))
		B <- n / (m - 1) *((mean(out[[1]]$calibration$s2_beta_0) - Psi_bar)^2 + (mean(out[[2]]$calibration$s2_beta_0) - Psi_bar)^2 + (mean(out[[3]]$calibration$s2_beta_0) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$calibration$s2_beta_0 - mean(out[[1]]$calibration$s2_beta_0))^2) + 1 / (n - 1) * sum((out[[2]]$calibration$s2_beta_0 - mean(out[[2]]$calibration$s2_beta_0))^2) + 1 / (n - 1) * sum((out[[3]]$calibration$s2_beta_0 - mean(out[[3]]$calibration$s2_beta_0))^2))
		R_hat_beta_0 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_beta_1
		Psi_bar <- mean(cbind(mean(out[[1]]$calibration$mu_beta_1), mean(out[[2]]$calibration$mu_beta_1), mean(out[[3]]$calibration$mu_beta_1)))
		B <- n / (m - 1) *((mean(out[[1]]$calibration$mu_beta_1) - Psi_bar)^2 + (mean(out[[2]]$calibration$mu_beta_1) - Psi_bar)^2 + (mean(out[[3]]$calibration$mu_beta_1) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$calibration$mu_beta_1 - mean(out[[1]]$calibration$mu_beta_1))^2) + 1 / (n - 1) * sum((out[[2]]$calibration$mu_beta_1 - mean(out[[2]]$calibration$mu_beta_1))^2) + 1 / (n - 1) * sum((out[[3]]$calibration$mu_beta_1 - mean(out[[3]]$calibration$mu_beta_1))^2))
		R_hat_beta_1 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_s2_beta_1
		Psi_bar <- mean(cbind(mean(out[[1]]$calibration$s2_beta_1), mean(out[[2]]$calibration$s2_beta_1), mean(out[[3]]$calibration$s2_beta_1)))
		B <- n / (m - 1) *((mean(out[[1]]$calibration$s2_beta_1) - Psi_bar)^2 + (mean(out[[2]]$calibration$s2_beta_1) - Psi_bar)^2 + (mean(out[[3]]$calibration$s2_beta_1) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$calibration$s2_beta_1 - mean(out[[1]]$calibration$s2_beta_1))^2) + 1 / (n - 1) * sum((out[[2]]$calibration$s2_beta_1 - mean(out[[2]]$calibration$s2_beta_1))^2) + 1 / (n - 1) * sum((out[[3]]$calibration$s2_beta_1 - mean(out[[3]]$calibration$s2_beta_1))^2))
		R_hat_beta_1 <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
	}
	
	## R_hat for W
	Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$W), rowMeans(out[[2]]$W), rowMeans(out[[3]]$W)))
	B <- n / (m - 1) *((rowMeans(out[[1]]$W) - Psi_bar)^2 + (rowMeans(out[[2]]$W) - Psi_bar)^2 + (rowMeans(out[[3]]$W) - Psi_bar)^2)
	W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$W - rowMeans(out[[1]]$W))^2) + 1 / (n - 1) * rowSums((out[[2]]$W - rowMeans(out[[2]]$W))^2) + 1 / (n - 1) * rowSums((out[[3]]$W - rowMeans(out[[3]]$W))^2))
	R_hat_W <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
	
	##
	## VS_lite model
	##
	if (model == 'VS'){
		
		## R_hat_Temp_min
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$Temp_min), rowMeans(out[[2]]$Temp_min), rowMeans(out[[3]]$Temp_min)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$Temp_min) - Psi_bar)^2 + (rowMeans(out[[2]]$Temp_min) - Psi_bar)^2 + (rowMeans(out[[3]]$Temp_min) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$Temp_min - rowMeans(out[[1]]$Temp_min))^2) + 1 / (n - 1) * rowSums((out[[2]]$Temp_min - rowMeans(out[[2]]$Temp_min))^2) + 1 / (n - 1) * rowSums((out[[3]]$Temp_min - rowMeans(out[[3]]$Temp_min))^2))
		R_hat_Temp_min <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_Temp_max
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$Temp_max), rowMeans(out[[2]]$Temp_max), rowMeans(out[[3]]$Temp_max)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$Temp_max) - Psi_bar)^2 + (rowMeans(out[[2]]$Temp_max) - Psi_bar)^2 + (rowMeans(out[[3]]$Temp_max) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$Temp_max - rowMeans(out[[1]]$Temp_max))^2) + 1 / (n - 1) * rowSums((out[[2]]$Temp_max - rowMeans(out[[2]]$Temp_max))^2) + 1 / (n - 1) * rowSums((out[[3]]$Temp_max - rowMeans(out[[3]]$Temp_max))^2))
		R_hat_Temp_max <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_P_min
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$P_min), rowMeans(out[[2]]$P_min), rowMeans(out[[3]]$P_min)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$P_min) - Psi_bar)^2 + (rowMeans(out[[2]]$P_min) - Psi_bar)^2 + (rowMeans(out[[3]]$P_min) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$P_min - rowMeans(out[[1]]$P_min))^2) + 1 / (n - 1) * rowSums((out[[2]]$P_min - rowMeans(out[[2]]$P_min))^2) + 1 / (n - 1) * rowSums((out[[3]]$P_min - rowMeans(out[[3]]$P_min))^2))
		R_hat_P_min <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_P_max
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$P_max), rowMeans(out[[2]]$P_max), rowMeans(out[[3]]$P_max)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$P_max) - Psi_bar)^2 + (rowMeans(out[[2]]$P_max) - Psi_bar)^2 + (rowMeans(out[[3]]$P_max) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$P_max - rowMeans(out[[1]]$P_max))^2) + 1 / (n - 1) * rowSums((out[[2]]$P_max - rowMeans(out[[2]]$P_max))^2) + 1 / (n - 1) * rowSums((out[[3]]$P_max - rowMeans(out[[3]]$P_max))^2))
		R_hat_P_max <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## output
		if(hierarchical){
			list(R_hat_beta_0 = R_hat_beta_0, R_hat_mu_beta_0 = R_hat_mu_beta_0, R_hat_s2_beta_0 = R_hat_s2_beta_0, R_hat_beta_1 = R_hat_beta_1, R_hat_mu_beta_1 = R_hat_mu_beta_1, R_hat_s2_beta_1 = R_hat_s2_beta_1, R_hat_s2 = R_hat_s2, R_hat_W = R_hat_W, R_hat_Temp_min = R_hat_Temp_min, R_hat_Temp_max = R_hat_Temp_max, R_hat_P_min = R_hat_P_min, R_hat_P_max = R_hat_P_max)	
		} else {
			list(R_hat_beta_0 = R_hat_beta_0, R_hat_beta_1 = R_hat_beta_1, R_hat_s2 = R_hat_s2, R_hat_W = R_hat_W, R_hat_Temp_min = R_hat_Temp_min, R_hat_Temp_max = R_hat_Temp_max, R_hat_P_min = R_hat_P_min, R_hat_P_max = R_hat_P_max)	
		}
	}
	
	##
	## Probit Model
	##
	
	else if (model == 'Probit'){
		## R_hat_gamma_T
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$gamma_T_save), rowMeans(out[[2]]$gamma_T_save), rowMeans(out[[3]]$gamma_T_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$gamma_T_save) - Psi_bar)^2 + (rowMeans(out[[2]]$gamma_T_save) - Psi_bar)^2 + (rowMeans(out[[3]]$gamma_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$gamma_T_save - rowMeans(out[[1]]$gamma_T_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$gamma_T_save - rowMeans(out[[2]]$gamma_T_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$gamma_T_save - rowMeans(out[[3]]$gamma_T_save))^2))
		R_hat_gamma_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_gamma_T
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_gamma_T_save), mean(out[[2]]$mu_gamma_T_save), mean(out[[3]]$mu_gamma_T_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_gamma_T_save) - Psi_bar)^2 + (mean(out[[2]]$mu_gamma_T_save) - Psi_bar)^2 + (mean(out[[3]]$mu_gamma_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_gamma_T_save - mean(out[[1]]$mu_gamma_T_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_gamma_T_save - mean(out[[2]]$mu_gamma_T_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_gamma_T_save - mean(out[[3]]$mu_gamma_T_save))^2))
		R_hat_mu_gamma_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_xi2_T
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$xi2_T_save), rowMeans(out[[2]]$xi2_T_save), rowMeans(out[[3]]$xi2_T_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$xi2_T_save) - Psi_bar)^2 + (rowMeans(out[[2]]$xi2_T_save) - Psi_bar)^2 + (rowMeans(out[[3]]$xi2_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$xi2_T_save - rowMeans(out[[1]]$xi2_T_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$xi2_T_save - rowMeans(out[[2]]$xi2_T_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$xi2_T_save - rowMeans(out[[3]]$xi2_T_save))^2))
		R_hat_xi2_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_xi2_T
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_xi2_T_save), mean(out[[2]]$mu_xi2_T_save), mean(out[[3]]$mu_xi2_T_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_xi2_T_save) - Psi_bar)^2 + (mean(out[[2]]$mu_xi2_T_save) - Psi_bar)^2 + (mean(out[[3]]$mu_xi2_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_xi2_T_save - mean(out[[1]]$mu_xi2_T_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_xi2_T_save - mean(out[[2]]$mu_xi2_T_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_xi2_T_save - mean(out[[3]]$mu_xi2_T_save))^2))
		R_hat_mu_xi2_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_gamma_P
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$gamma_P_save), rowMeans(out[[2]]$gamma_P_save), rowMeans(out[[3]]$gamma_P_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$gamma_P_save) - Psi_bar)^2 + (rowMeans(out[[2]]$gamma_P_save) - Psi_bar)^2 + (rowMeans(out[[3]]$gamma_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$gamma_P_save - rowMeans(out[[1]]$gamma_P_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$gamma_P_save - rowMeans(out[[2]]$gamma_P_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$gamma_P_save - rowMeans(out[[3]]$gamma_P_save))^2))
		R_hat_gamma_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_gamma_P
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_gamma_P_save), mean(out[[2]]$mu_gamma_P_save), mean(out[[3]]$mu_gamma_P_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_gamma_P_save) - Psi_bar)^2 + (mean(out[[2]]$mu_gamma_P_save) - Psi_bar)^2 + (mean(out[[3]]$mu_gamma_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_gamma_P_save - mean(out[[1]]$mu_gamma_P_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_gamma_P_save - mean(out[[2]]$mu_gamma_P_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_gamma_P_save - mean(out[[3]]$mu_gamma_P_save))^2))
		R_hat_mu_gamma_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_xi2_P
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$xi2_P_save), rowMeans(out[[2]]$xi2_P_save), rowMeans(out[[3]]$xi2_P_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$xi2_P_save) - Psi_bar)^2 + (rowMeans(out[[2]]$xi2_P_save) - Psi_bar)^2 + (rowMeans(out[[3]]$xi2_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$xi2_P_save - rowMeans(out[[1]]$xi2_P_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$xi2_P_save - rowMeans(out[[2]]$xi2_P_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$xi2_P_save - rowMeans(out[[3]]$xi2_P_save))^2))
		R_hat_xi2_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_xi2_P
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_xi2_P_save), mean(out[[2]]$mu_xi2_P_save), mean(out[[3]]$mu_xi2_P_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_xi2_P_save) - Psi_bar)^2 + (mean(out[[2]]$mu_xi2_P_save) - Psi_bar)^2 + (mean(out[[3]]$mu_xi2_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_xi2_P_save - mean(out[[1]]$mu_xi2_P_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_xi2_P_save - mean(out[[2]]$mu_xi2_P_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_xi2_P_save - mean(out[[3]]$mu_xi2_P_save))^2))
		R_hat_mu_xi2_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## write output
		if(hierarchical){
			list(R_hat_beta_0 = R_hat_beta_0, R_hat_mu_beta_0 = R_hat_mu_beta_0, R_hat_s2_beta_0 = R_hat_s2_beta_0, R_hat_beta_1 = R_hat_beta_1, R_hat_mu_beta_1 = R_hat_mu_beta_1, R_hat_s2_beta_1 = R_hat_s2_beta_1, R_hat_s2 = R_hat_s2, R_hat_W = R_hat_W, R_hat_gamma_T = R_hat_gamma_T, R_hat_xi2_T = R_hat_xi2_T, R_hat_gamma_P = R_hat_gamma_P, R_hat_xi2_P = R_hat_xi2_P, R_hat_mu_gamma_T = R_hat_mu_gamma_T, R_hat_mu_xi2_T = R_hat_mu_xi2_T, R_hat_mu_gamma_P = R_hat_mu_gamma_P, R_hat_mu_xi2_P = R_hat_mu_xi2_P)	
		} else {
		list(R_hat_beta_0 = R_hat_beta_0, R_hat_beta_1 = R_hat_beta_1, R_hat_s2 = R_hat_s2, R_hat_W = R_hat_W, R_hat_gamma_T = R_hat_gamma_T, R_hat_xi2_T = R_hat_xi2_T, R_hat_gamma_P = R_hat_gamma_P, R_hat_xi2_P = R_hat_xi2_P, R_hat_mu_gamma_T = R_hat_mu_gamma_T, R_hat_mu_xi2_T = R_hat_mu_xi2_T, R_hat_mu_gamma_P = R_hat_mu_gamma_P, R_hat_mu_xi2_P = R_hat_mu_xi2_P)	
		}
	}
	
	##
	## Mixed model
	##
	
	else if (model == 'Mixed'){
		
		## R_hat_Temp_min
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$Temp_min_save), rowMeans(out[[2]]$Temp_min_save), rowMeans(out[[3]]$Temp_min_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$Temp_min_save) - Psi_bar)^2 + (rowMeans(out[[2]]$Temp_min_save) - Psi_bar)^2 + (rowMeans(out[[3]]$Temp_min_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$Temp_min_save - rowMeans(out[[1]]$Temp_min_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$Temp_min_save - rowMeans(out[[2]]$Temp_min_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$Temp_min_save - rowMeans(out[[3]]$Temp_min_save))^2))
		R_hat_Temp_min <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_Temp_max
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$Temp_max_save), rowMeans(out[[2]]$Temp_max_save), rowMeans(out[[3]]$Temp_max_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$Temp_max_save) - Psi_bar)^2 + (rowMeans(out[[2]]$Temp_max_save) - Psi_bar)^2 + (rowMeans(out[[3]]$Temp_max_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$Temp_max_save - rowMeans(out[[1]]$Temp_max_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$Temp_max_save - rowMeans(out[[2]]$Temp_max_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$Temp_max_save - rowMeans(out[[3]]$Temp_max_save))^2))
		R_hat_Temp_max <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_P_min
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$P_min_save), rowMeans(out[[2]]$P_min_save), rowMeans(out[[3]]$P_min_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$P_min_save) - Psi_bar)^2 + (rowMeans(out[[2]]$P_min_save) - Psi_bar)^2 + (rowMeans(out[[3]]$P_min_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$P_min_save - rowMeans(out[[1]]$P_min_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$P_min_save - rowMeans(out[[2]]$P_min_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$P_min_save - rowMeans(out[[3]]$P_min_save))^2))
		R_hat_P_min <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_P_max
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$P_max_save), rowMeans(out[[2]]$P_max_save), rowMeans(out[[3]]$P_max_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$P_max_save) - Psi_bar)^2 + (rowMeans(out[[2]]$P_max_save) - Psi_bar)^2 + (rowMeans(out[[3]]$P_max_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$P_max_save - rowMeans(out[[1]]$P_max_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$P_max_save - rowMeans(out[[2]]$P_max_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$P_max_save - rowMeans(out[[3]]$P_max_save))^2))
		R_hat_P_max <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_gamma_T
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$gamma_T_save), rowMeans(out[[2]]$gamma_T_save), rowMeans(out[[3]]$gamma_T_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$gamma_T_save) - Psi_bar)^2 + (rowMeans(out[[2]]$gamma_T_save) - Psi_bar)^2 + (rowMeans(out[[3]]$gamma_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$gamma_T_save - rowMeans(out[[1]]$gamma_T_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$gamma_T_save - rowMeans(out[[2]]$gamma_T_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$gamma_T_save - rowMeans(out[[3]]$gamma_T_save))^2))
		R_hat_gamma_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_gamma_T
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_gamma_T_save), mean(out[[2]]$mu_gamma_T_save), mean(out[[3]]$mu_gamma_T_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_gamma_T_save) - Psi_bar)^2 + (mean(out[[2]]$mu_gamma_T_save) - Psi_bar)^2 + (mean(out[[3]]$mu_gamma_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_gamma_T_save - mean(out[[1]]$mu_gamma_T_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_gamma_T_save - mean(out[[2]]$mu_gamma_T_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_gamma_T_save - mean(out[[3]]$mu_gamma_T_save))^2))
		R_hat_mu_gamma_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_xi2_T
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$xi2_T_save), rowMeans(out[[2]]$xi2_T_save), rowMeans(out[[3]]$xi2_T_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$xi2_T_save) - Psi_bar)^2 + (rowMeans(out[[2]]$xi2_T_save) - Psi_bar)^2 + (rowMeans(out[[3]]$xi2_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$xi2_T_save - rowMeans(out[[1]]$xi2_T_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$xi2_T_save - rowMeans(out[[2]]$xi2_T_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$xi2_T_save - rowMeans(out[[3]]$xi2_T_save))^2))
		R_hat_xi2_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_xi2_T
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_xi2_T_save), mean(out[[2]]$mu_xi2_T_save), mean(out[[3]]$mu_xi2_T_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_xi2_T_save) - Psi_bar)^2 + (mean(out[[2]]$mu_xi2_T_save) - Psi_bar)^2 + (mean(out[[3]]$mu_xi2_T_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_xi2_T_save - mean(out[[1]]$mu_xi2_T_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_xi2_T_save - mean(out[[2]]$mu_xi2_T_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_xi2_T_save - mean(out[[3]]$mu_xi2_T_save))^2))
		R_hat_mu_xi2_T <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_gamma_P
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$gamma_P_save), rowMeans(out[[2]]$gamma_P_save), rowMeans(out[[3]]$gamma_P_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$gamma_P_save) - Psi_bar)^2 + (rowMeans(out[[2]]$gamma_P_save) - Psi_bar)^2 + (rowMeans(out[[3]]$gamma_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$gamma_P_save - rowMeans(out[[1]]$gamma_P_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$gamma_P_save - rowMeans(out[[2]]$gamma_P_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$gamma_P_save - rowMeans(out[[3]]$gamma_P_save))^2))
		R_hat_gamma_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_gamma_P
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_gamma_P_save), mean(out[[2]]$mu_gamma_P_save), mean(out[[3]]$mu_gamma_P_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_gamma_P_save) - Psi_bar)^2 + (mean(out[[2]]$mu_gamma_P_save) - Psi_bar)^2 + (mean(out[[3]]$mu_gamma_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_gamma_P_save - mean(out[[1]]$mu_gamma_P_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_gamma_P_save - mean(out[[2]]$mu_gamma_P_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_gamma_P_save - mean(out[[3]]$mu_gamma_P_save))^2))
		R_hat_mu_gamma_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_xi2_P
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$xi2_P_save), rowMeans(out[[2]]$xi2_P_save), rowMeans(out[[3]]$xi2_P_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$xi2_P_save) - Psi_bar)^2 + (rowMeans(out[[2]]$xi2_P_save) - Psi_bar)^2 + (rowMeans(out[[3]]$xi2_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$xi2_P_save - rowMeans(out[[1]]$xi2_P_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$xi2_P_save - rowMeans(out[[2]]$xi2_P_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$xi2_P_save - rowMeans(out[[3]]$xi2_P_save))^2))
		R_hat_xi2_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_mu_xi2_P
		Psi_bar <- mean(cbind(mean(out[[1]]$mu_xi2_P_save), mean(out[[2]]$mu_xi2_P_save), mean(out[[3]]$mu_xi2_P_save)))
		B <- n / (m - 1) *((mean(out[[1]]$mu_xi2_P_save) - Psi_bar)^2 + (mean(out[[2]]$mu_xi2_P_save) - Psi_bar)^2 + (mean(out[[3]]$mu_xi2_P_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * sum((out[[1]]$mu_xi2_P_save - mean(out[[1]]$mu_xi2_P_save))^2) + 1 / (n - 1) * sum((out[[2]]$mu_xi2_P_save - mean(out[[2]]$mu_xi2_P_save))^2) + 1 / (n - 1) * sum((out[[3]]$mu_xi2_P_save - mean(out[[3]]$mu_xi2_P_save))^2))
		R_hat_mu_xi2_P <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## R_hat_x
		Psi_bar <- rowMeans(cbind(rowMeans(out[[1]]$x_save), rowMeans(out[[2]]$x_save), rowMeans(out[[3]]$x_save)))
		B <- n / (m - 1) *((rowMeans(out[[1]]$x_save) - Psi_bar)^2 + (rowMeans(out[[2]]$x_save) - Psi_bar)^2 + (rowMeans(out[[3]]$x_save) - Psi_bar)^2)
		W <- 1 / m * (1 / (n - 1) * rowSums((out[[1]]$x_save - rowMeans(out[[1]]$x_save))^2) + 1 / (n - 1) * rowSums((out[[2]]$x_save - rowMeans(out[[2]]$x_save))^2) + 1 / (n - 1) * rowSums((out[[3]]$x_save - rowMeans(out[[3]]$x_save))^2))
		R_hat_x <- sqrt(((n - 1) / n * W + 1 / n * B) / W)
		
		## write output
		list(R_hat_beta_0 = R_hat_beta_0, R_hat_beta_1 = R_hat_beta_1,  R_hat_s2 = R_hat_s2, R_hat_W = R_hat_W, R_hat_Temp_min = R_hat_Temp_min, R_hat_Temp_max = R_hat_Temp_max, R_hat_P_min = R_hat_P_min, R_hat_P_max = R_hat_P_max, R_hat_gamma_T = R_hat_gamma_T, R_hat_xi2_T = R_hat_xi2_T, R_hat_gamma_P = R_hat_gamma_P, R_hat_xi2_P = R_hat_xi2_P, R_hat_mu_gamma_T = R_hat_mu_gamma_T, R_hat_mu_xi2_T = R_hat_mu_xi2_T, R_hat_mu_gamma_P = R_hat_mu_gamma_P, R_hat_mu_xi2_P = R_hat_mu_xi2_P, R_hat_x = R_hat_x)
	}
}