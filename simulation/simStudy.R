set.seed(101)

##
## Libraries and functions
##

library(myFunctions)
library(bayesTreeRing)

i = 3
##
## select simulation model
##

if (i == 1){ ## sim using VS
	model = 'VS'
} else if (i == 2){ ## sim using Probit
	model = 'Probit'
} else if(i == 3){ ## sim using Mixed
	model = 'Mixed'
}


##
## select trend
## 

trend = TRUE

##
## s2_noise
##

s2_noise <- 0.5

##
## Simulate data
##

dataSim <- simData(model, s2_noise, trend)

##
## Fit empirical Bayesian climate process model	
## 

process <- makeMCARTrend(dataSim$W[1:12, 447:556], dataSim$W[13:24, 447:556])

n_mcmc <- 15000

##
## select model fit
##

##
## VS Lite Model
##

cat('\n Iteration',  i, 'VS model fit \n \n')
source('~/Mechanistic-Tree-Ring/simulation/fit/fitVS.R')
makeVSHierarchical(dataSim, model, process, trend, n_mcmc)

##
## fit Probit model
##
cat('\n Iteration',  i, 'Probit model fit \n \n')
source('~/Mechanistic-Tree-Ring/simulation/fit/fitProbit.R')
makeProbitHierarchicalRandom(dataSim, model, process, trend, n_mcmc)

##
## fit Mixed model
##
cat('\n Iteration',  i, 'Mixed model fit \n \n')
source('~/Mechanistic-Tree-Ring/simulation/fit/fitMixed.R')
makeMixedHierarchicalRandom(dataSim, model, process, trend, n_mcmc)

