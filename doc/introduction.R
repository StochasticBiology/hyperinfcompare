## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(hyperinf)
library(hyperinfcompare)
library(ggpubr)

## -----------------------------------------------------------------------------
# create a very simple EvAM dataset
set.seed(1)
n = 4
data.a = matrix(rep(c(0,0,0,0,1, 
                      0,0,0,1,1, 
                      0,0,1,1,1,
                      0,1,1,1,1), n/4), byrow = TRUE, ncol=5, nrow=n)

ggarrange(plot_hyperinf_data(data.a), 
          plot_hyperinf_data(1-data.a))

## -----------------------------------------------------------------------------
# inference using the dataset and its "inverse" -- by default using HyperHMM
fit.1a = hyperinf(data.a, nboot = 30)
fit.2a = hyperinf(1-data.a, nboot = 30)

# plot the comparison
plot.co.a = plot_hyperinf_compare_orderings(fit.1a, fit.2a, 
                                            expt.names = c("i", "ii"),
                                            thetastep=3, p.scale=0.3)
print(plot.co.a)

# the comparison itself
compare_orderings(fit.1a, fit.2a)

