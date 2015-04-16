
##
## Libraries and Functions
##

library(xtable)

sumCRPS <- function(out){
	c(sum(out$crps_T), sum(out$crps_T_summer), sum(out$crps_P), sum(out$crps_P_summer))
}

sumCRPSClimate <- function(out){
	c(sum(out$climate_CRPS_T), sum(out$climate_CRPS_T_summer), sum(out$climate_CRPS_P), sum(out$climate_CRPS_P_summer))
}


meanCRPS <- function(out){
	c(mean(out$crps_T), mean(out$crps_T_summer), mean(out$crps_P), mean(out$crps_P_summer))
}

meanCRPSClimate <- function(out){
	c(mean(out$climate_CRPS_T), mean(out$climate_CRPS_T_summer), mean(out$climate_CRPS_P), mean(out$climate_CRPS_P_summer))
}


histPIT <- function(out){
	layout(matrix(1:4, 2, 2))
	hist(out$PIT_T)
	hist(out$PIT_T_summer)
	hist(out$PIT_P)
	hist(out$PIT_P_summer)
}


analyzeResults <- function(i){
	
	source('~/Mechanistic-Tree-Ring/simulation/makeRHatFull.R')
	source('~/Mechanistic-Tree-Ring/simulation/evaluatePredictions.R')
	
	##
	## Data simulated with VS-Lite model
	##
	
	# 	R_hat <- vector("list")
	# 	eval_sim <- vector("list")
	# 	simRun <- vector("list")
	
	## Hierarchical model fit
	## VS-Lite Model fit
	if(i == 1){
		load('~/Mechanistic-Tree-Ring/data/simVSfitVSHierarchical.RData')
		R_hat <- make_R_hat(simVSfitVSHierarchical, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simVSfitVSHierarchical, dataSim, model = 'VS')
		simRun <- "simVSfitVSHierarchical"
	}
	## Hierarchical model fit
	## VS-Lite Mis-specified Model fit
	if(i == 2){
		load('~/Mechanistic-Tree-Ring/data/simVSfitVSHierarchicalMisSpec.RData')
		R_hat <- make_R_hat(simVSfitVSHierarchicalMisSpec, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simVSfitVSHierarchicalMisSpec, dataSim, model = 'VS')
		simRun <- "simVSfitVSHierarchicalMisSpec"
	}
	## Hierarchical Random model fit
	## Probit Model fit
	if(i == 3){
		load('~/Mechanistic-Tree-Ring/data/simVSfitProbitHierarchicalRandom.RData')
		R_hat <- make_R_hat(simVSfitProbitHierarchicalRandom, model = 'Probit', hierarchical = TRUE, random = TRUE)
		# R_hat
		eval <- evaluatePredictions(simVSfitProbitHierarchicalRandom, dataSim, model = 'Probit')
		simRun <- "simVSfitProbitHierarchicalRandom"
	}
	if(i == 4){	
		## Mixed Model fit
		load('~/Mechanistic-Tree-Ring/data/simVSfitMixedHierarchicalRandom.RData')
		R_hat <- make_R_hat(simVSfitMixedHierarchicalRandom, model = 'Mixed', hierarchical = TRUE, random = TRUE)
		eval <- evaluatePredictions(simVSfitMixedHierarchicalRandom, dataSim, model = 'Mixed')
		simRun <- "simVSfitMixedHierarchicalRandom"
	}
	##
	## Data simulated with Probit model
	##
	
	## Hierarchical model fit
	## VS-Lite Model fit
	if(i == 5){
		load('~/Mechanistic-Tree-Ring/data/simProbitfitVSHierarchical.RData')
		R_hat <- make_R_hat(simProbitfitVSHierarchical, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simProbitfitVSHierarchical, dataSim, model = 'VS')
		simRun <- "simProbitfitVSHierarchical"
	}

	## Hierarchical model fit
	## VS-Lite Mis-specified Model fit
	if(i == 6){
		load('~/Mechanistic-Tree-Ring/data/simProbitfitVSHierarchicalMisSpec.RData')
		R_hat <- make_R_hat(simProbitfitVSHierarchicalMisSpec, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simProbitfitVSHierarchicalMisSpec, dataSim, model = 'VS')
		simRun <- "simProbitfitVSHierarchicalMisSpec"
	}
	
	## Hierarchical Random model fit
	## Probit Model fit
	if(i == 7){
		load('~/Mechanistic-Tree-Ring/data/simProbitfitProbitHierarchicalRandom.RData')
		R_hat <- make_R_hat(simProbitfitProbitHierarchicalRandom, model = 'Probit', hierarchical = TRUE, random = TRUE)
		# R_hat
		eval <- evaluatePredictions(simProbitfitProbitHierarchicalRandom, dataSim, model = 'Probit')
		simRun <- "simProbitfitProbitHierarchicalRandom"
	}
	if(i == 8){	
		## Mixed Model fit
		load('~/Mechanistic-Tree-Ring/data/simProbitfitMixedHierarchicalRandom.RData')
		R_hat <- make_R_hat(simProbitfitMixedHierarchicalRandom, model = 'Mixed', hierarchical = TRUE, random = TRUE)
		eval <- evaluatePredictions(simProbitfitMixedHierarchicalRandom, dataSim, model = 'Mixed')
		simRun <- "simProbitfitMixedHierarchicalRandom"
	}
	## 
	## Data simulated with Mixed model
	##
	
	## Hierarchical model fit
	## VS-Lite Model fit
	if(i == 9){
		load('~/Mechanistic-Tree-Ring/data/simMixedfitVSHierarchical.RData')
		R_hat <- make_R_hat(simMixedfitVSHierarchical, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simMixedfitVSHierarchical, dataSim, model = 'VS')
		simRun <- "simMixedfitVSHierarchical"
	}
	
	
	## Hierarchical model fit
	## VS-Lite Mis-specified Model fit
	if(i == 10){
		load('~/Mechanistic-Tree-Ring/data/simMixedfitVSHierarchicalMisSpec.RData')
		R_hat <- make_R_hat(simMixedfitVSHierarchicalMisSpec, model = 'VS', hierarchical = TRUE)
		eval <- evaluatePredictions(simMixedfitVSHierarchicalMisSpec, dataSim, model = 'VS')
		simRun <- "simMixedfitVSHierarchicalMisSpec"
	}
	
	## Hierarchical Random model fit
	## Probit Model fit
	if(i == 11){
		load('~/Mechanistic-Tree-Ring/data/simMixedfitProbitHierarchicalRandom.RData')
		R_hat <- make_R_hat(simMixedfitProbitHierarchicalRandom, model = 'Probit', hierarchical = TRUE, random = TRUE)
		# R_hat
		eval <- evaluatePredictions(simMixedfitProbitHierarchicalRandom, dataSim, model = 'Probit')
		simRun <- "simMixedfitProbitHierarchicalRandom"
	}
	if(i == 12){	
		## Mixed Model fit
		load('~/Mechanistic-Tree-Ring/data/simMixedfitMixedHierarchicalRandom.RData')
		R_hat <- make_R_hat(simMixedfitMixedHierarchicalRandom, model = 'Mixed', hierarchical = TRUE, random = TRUE)
		eval <- evaluatePredictions(simMixedfitMixedHierarchicalRandom, dataSim, model = 'Mixed')
		simRun <- "simMixedfitMixedHierarchicalRandom"
	}
	return(list(R_hat = R_hat, eval = eval, simRun = simRun))
}

out_analysis <- lapply(1:12, analyzeResults)

# save.image('~/Mechanistic-Tree-Ring/data/simResultsMar11.RData')
# load('~/Mechanistic-Tree-Ring/data/simResultsMar3.RData')

checkR_hat <- function(i){
	max(unlist(out_analysis[[i]]$R_hat))
}

sapply(1:12, checkR_hat)

wrapResults <- function(i){
		if(i %in% c(1, 5, 9)){
			cbind(sumCRPSClimate(out_analysis[[i]]$eval),
						sumCRPS(out_analysis[[i]]$eval))
# 			cbind(meanCRPSClimate(out_analysis[[i]]$eval),
# 						meanCRPS(out_analysis[[i]]$eval))
		} else {
			sumCRPS(out_analysis[[i]]$eval)
# 			meanCRPS(out_analysis[[i]]$eval)
		}
}

resultsVSSim <- data.frame(lapply(1:4, wrapResults))
rownames(resultsVSSim) <- c('crps_T', 'crps_T_summer', 'crps_P', 'crps_P_summer')
colnames(resultsVSSim) <- c('Clim', 'VS', 'VS-Mis', 'P', 'M')


resultsPSim <- data.frame(lapply(5:8, wrapResults))
rownames(resultsPSim) <- c('crps_T', 'crps_T_summer', 'crps_P', 'crps_P_summer')
colnames(resultsPSim) <- c('Clim', 'VS', 'VS-Mis', 'P', 'M')


resultsMSim <- data.frame(lapply(9:12, wrapResults))
rownames(resultsMSim) <- c('crps_T', 'crps_T_summer', 'crps_P', 'crps_P_summer')
colnames(resultsMSim) <- c('Clim', 'VS', 'VS-Mis', 'P', 'M')



xtable(resultsVSSim, digits = 2)
xtable(resultsPSim, digits = 2)
xtable(resultsMSim, digits = 2)
##
## results Table
##

## Sim model    |            VS            |            P            |           N/P           | ...
##              |                          |                         |                         | ...
# Fit model     |                          |                         |                         | ...
##              | Clim | VS | Pr | N/P | M | Clim | VS | P | N/P | M | Clim | VS | P | N/P | M | ...      
## crps_T       |                          |                         |                         | ...
## crps_P       |                          |                         |                         | ...
## crps_T_summer|                          |                         |                         | ...
## crps_P_summer|                          |                         |                         | ...



