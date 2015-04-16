evaluatePredictions <- function(out, data, model, n_thin = 5){
	
  library(myFunctions)
  library(bayesTreeRing)
  # 	cat("plotting predictions \n")
    n_iter <- dim(out[[1]]$calibration$beta_0)[2]
    thin_idx <- seq(1, n_iter, by = n_thin)
    year <- 2005 - (556:1)
	
    W_T_month <- cbind(downscaleMonth(out[[1]]$W[1:12, , thin_idx]),
                       downscaleMonth(out[[2]]$W[1:12, , thin_idx]),
                       downscaleMonth(out[[3]]$W[1:12, , thin_idx]))
    W_T_summer <- cbind(downscaleMonth(out[[1]]$W[4:9, , thin_idx]), 
                        downscaleMonth(out[[2]]$W[4:9, , thin_idx]), 
                        downscaleMonth(out[[3]]$W[4:9, , thin_idx]))
    W_T_true <- apply(data$W[1:12, ], 2, mean)
    W_T_true_summer <- apply(data$W[4:9, ], 2, mean)
    W_P_month <- cbind(downscaleMonth(out[[1]]$W[13:24, , thin_idx]),
                       downscaleMonth(out[[2]]$W[13:24, , thin_idx]),
                       downscaleMonth(out[[3]]$W[13:24, , thin_idx]))
    W_P_summer <- cbind(downscaleMonth(out[[1]]$W[16:21, , thin_idx]), 
                        downscaleMonth(out[[2]]$W[16:21, , thin_idx]), 
                        downscaleMonth(out[[3]]$W[16:21, , thin_idx]))
    W_P_true <- apply(data$W[13:24, ], 2, mean)
    W_P_true_summer <- apply(data$W[16:21, ], 2, mean)
    
	##
	## Annual climate reconstructions
	##
    
# 	layout(matrix(1:2, 2, 1))
# 	## Temperature precipitation plots
# 	matplot(apply(W_T_month, 1, median), type = 'l', ylim = c(6, 12),
# 					main = paste(model, 'Temp'), ylab = 'T', xaxt = 'n', xlab = 'Year')
# 	axis(1, at = 1:length(year), labels = year)
# 	matplot(apply(W_T_month, 1, quantile, prob = 0.975), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(apply(W_T_month, 1, quantile, prob = 0.025), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(W_T_true, type = 'l', add = TRUE, col = 'blue')
#   ## Log precipitation plots
# 	matplot(apply(W_P_month, 1, median), type = 'l', ylim = c(4, 5),
# 					main = paste(model, 'Log Precip'), ylab = 'log(P)', xaxt = 'n', 
# 					xlab = 'Year')
# 	axis(1, at = 1:length(year), labels = year)
# 	matplot(apply(W_P_month, 1, quantile, prob = 0.975), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(apply(W_P_month, 1, quantile, prob = 0.025), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(W_P_true, type = 'l', add = TRUE, col = 'blue')
# 	
# 	##
# 	## Summer plots
# 	##
# 	## Summer Temperature plots
# 	matplot(apply(W_T_summer, 1, median), type = 'l', ylim = c(14, 19),
# 					main = paste(model, 'Summer Temp'), ylab = 'T', xaxt = 'n', xlab = 'Year')
# 	axis(1, at = 1:length(year), labels = year)
# 	matplot(apply(W_T_summer, 1, quantile, prob = 0.975), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(apply(W_T_summer, 1, quantile, prob = 0.025), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(W_T_true_summer, type = 'l', add = TRUE, col = 'blue')
# 	## Summer precipitation plots
# 	matplot(apply(W_P_summer, 1, median), type = 'l', ylim = c(4, 5),
# 					main = paste(model, 'Summer Log Precipitation'), ylab = 'T', xaxt = 'n', xlab = 'Year')
# 	axis(1, at = 1:length(year), labels = year)
# 	matplot(apply(W_P_summer, 1, quantile, prob = 0.975), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(apply(W_P_summer, 1, quantile, prob = 0.025), type = 'l', 
# 					col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
# 	matplot(W_P_true_summer, type = 'l', add = TRUE, col = 'blue')
	
    N_T <- 110
    t <- length(W_T_true)
    
    ##
    ## CRPS scores
    ##
    
    cat("CRPS scores for model \n")
    crps_T <- makeCRPS(t - N_T, W_T_month[ -1, ], W_T_true[ -1 ])
    crps_T_summer <- makeCRPS(t - N_T, W_T_summer[ -1, ], W_T_true_summer[ -1 ])
    crps_P <- makeCRPS(t - N_T, W_P_month[ - 1, ], W_P_true[ -1 ])
    crps_P_summer <- makeCRPS(t - N_T, W_P_summer[ -1, ], W_P_true_summer[ -1 ])
    
    if(model == "VS"){
        cat("CRPS scores for climatology \n")
        n_iter <- dim(out[[1]]$calibration$beta_0)[2] * 3 / n_thin
        Temp_hat_summer <- predictClimatology(n_iter, data$W[4:9, 447:556])
        P_hat_summer <- predictClimatology(n_iter, data$W[16:21, 447:556])
        Temp_hat <- predictClimatology(n_iter, data$W[1:12, 447:556])
        P_hat <- predictClimatology(n_iter, data$W[13:24, 447:556])
        climate_CRPS_T <- makeCRPS(445, Temp_hat[ - 1, ], W_T_true[ -1 ])	
        climate_CRPS_P <- makeCRPS(445, P_hat[ - 1, ], W_P_true[ -1 ])
        climate_CRPS_T_summer <- makeCRPS(445, Temp_hat_summer[ - 1, ], W_T_true_summer[ -1 ])	
        climate_CRPS_P_summer <- makeCRPS(445, P_hat_summer[ - 1, ], W_P_true_summer[ -1 ])
    }
    
    ##
    ## Probability Integral Transform
    ##
    
    cat("making PIT transforms \n")
    PIT_T <- makePIT(t - N_T, W_T_month, W_T_true)
    PIT_T_summer <- makePIT(t - N_T, W_T_summer, W_T_true_summer)
    PIT_P <- makePIT(t - N_T, W_P_month, W_P_true)
    PIT_P_summer <- makePIT(t - N_T, W_P_summer, W_P_true_summer)
    
    ##
    ## return values
    ##
    if(model == "VS"){
        list(crps_T = crps_T, crps_T_summer = crps_T_summer, crps_P = crps_P, crps_P_summer = crps_P_summer, climate_CRPS_T = climate_CRPS_T, climate_CRPS_P = climate_CRPS_P, climate_CRPS_T_summer = climate_CRPS_T_summer, climate_CRPS_P_summer = climate_CRPS_P_summer, PIT_T = PIT_T, PIT_T_summer = PIT_T_summer, PIT_P = PIT_P, PIT_P_summer = PIT_P_summer)
    } else {
        list(crps_T = crps_T, crps_T_summer = crps_T_summer, crps_P = crps_P, crps_P_summer = crps_P_summer, PIT_T = PIT_T, PIT_T_summer = PIT_T_summer, PIT_P = PIT_P, PIT_P_summer = PIT_P_summer)
    }		
}
