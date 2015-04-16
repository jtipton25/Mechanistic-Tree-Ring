
##
## Libraries and Functions
##

library(xtable)
library(bayesTreeRing)
	
source('~/Mechanistic-Tree-Ring/simulation/makeRHatFull.R')
source('~/Mechanistic-Tree-Ring/simulation/evaluatePredictions.R')
	
load('~/Mechanistic-Tree-Ring/data/fitMixedHierarchicalRandomHudson.RData')
R_hat <- make_R_hat(out, model = 'Mixed', hierarchical = TRUE)
max(unlist(R_hat))

## Define index for year labels
year <- 2005 - (556:1)
t <- 556-1

TempHudson <- cbind(downscaleMonth(out[[1]]$W[1:12, , ]),
										downscaleMonth(out[[2]]$W[1:12, , ]),
										downscaleMonth(out[[3]]$W[1:12, , ]))
PHudson <- cbind(downscaleMonth(out[[1]]$W[13:24, , ]),
								 downscaleMonth(out[[2]]$W[13:24, , ]),
								 downscaleMonth(out[[3]]$W[13:24, , ]))
# 	cbind(out[[1]]$W, out[[2]]$W, out[[3]]$W)



# tiff(file = '~/treeRing/plots/MixedHudsonReconstruction.tiff', width = 12, 
# 		 height = 6, res  = 600, units = 'in')
# jpeg(file = '~/treeRing/plots/MixedHudsonReconstruction.jpg', width = 12, 
# 		 height = 6, res  = 100, units = 'in')
layout(matrix(1:2, 2, 1))
matplot(apply(TempHudson, 1, median)[-1], type = 'n', ylim = c(6, 12),
				main = paste('Temperature reconstruction'), 
				ylab = expression(~degree~C), xaxt = 'n', xlab = 'Year')
matplot(c(rep(NA, 445), apply(TempHudson, 1, median)[-1][445:555]), type = 'l', add = TRUE)
axis(1, at = seq(2, 556, 50), labels = year[seq(2, 556, 50)])
for(i in 1:19){
	polygon(c(1:t, t:1), c(apply(TempHudson, 1, quantile, prob = (0.5 + i / 40))[-1],
												 rev(apply(TempHudson, 1, quantile, prob = (0.5 + (i - 1) / 40))[-1])), 
					col = adjustcolor('grey60', alpha.f = (1 - (i - 1) / 20)), border = NA)  
	polygon(c(1:t, t:1), c(apply(TempHudson, 1, quantile, prob = (0.5 - i / 40))[-1],
												 rev(apply(TempHudson, 1, quantile, prob = (0.5 - (i - 1) / 40))[-1])), 
					col = adjustcolor('grey60', alpha.f = (1 - (i - 1) / 20)), border = NA)  
}
matplot(apply(TempHudson, 1, quantile, prob = 0.975)[-1], type = 'l', 
				lty = 2, lwd = 0.25, add = TRUE)
matplot(apply(TempHudson, 1, quantile, prob = 0.025)[-1], type = 'l', 
				lty = 2, lwd = 0.25, add = TRUE)





matplot(apply(PHudson, 1, median)[-1], type = 'n', ylim = c(3.9, 4.9),
				main = paste('Log Precipitation reconstruction'), 
				ylab = expression(~degree~C), xaxt = 'n', xlab = 'Year')
axis(1, at = seq(2, 556, 50), labels = year[seq(2, 556, 50)])
matplot(c(rep(NA, 445), apply(PHudson, 1, median)[-1][445:555]), type = 'l', add = TRUE)
for(i in 1:19){
	polygon(c(1:t, t:1), c(apply(PHudson, 1, quantile, prob = (0.5 + i / 40))[-1],
												 rev(apply(PHudson, 1, quantile, prob = (0.5 + (i - 1) / 40))[-1])), 
					col = adjustcolor('grey60', alpha.f = (1 - (i - 1) / 20)), border = NA)  
	polygon(c(1:t, t:1), c(apply(PHudson, 1, quantile, prob = (0.5 - i / 40))[-1],
												 rev(apply(PHudson, 1, quantile, prob = (0.5 - (i - 1) / 40))[-1])), 
					col = adjustcolor('grey60', alpha.f = (1 - (i - 1) / 20)), border = NA)  
}
matplot(apply(PHudson, 1, quantile, prob = 0.975)[-1], type = 'l', 
				lty = 2, lwd = 0.25, add = TRUE)
matplot(apply(PHudson, 1, quantile, prob = 0.025)[-1], type = 'l', 
				lty = 2, lwd = 0.25, add = TRUE)

pederson <- read.table(file = '~/Google Drive/NYCDroughtPederson.txt',
											 sep = '', header = TRUE)
tmp <- rbind(matrix(NA, 82, 2, dimnames = list(c(), c("YEAR", "drought"))), 
						 pederson[1:474, ])
matplot((tmp$drought[-1] / sd(tmp$drought[!is.na(tmp$drought)])) * sd(PHudson[-1, ]) + 
					mean(PHudson[-1, ]), type = 'l', col = 'black', lty = 1, add = TRUE)
dev.off()
