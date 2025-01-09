########################################
#  Simulations to compare test results from gplsim package and gamFactory
########################################
rm(list=ls())
library(gplsim)
library(gamFactory)
library(parallel)
library(ggplot2)
library(tidyverse)
source("./simulations_gplsim_functions.R")

# parameter settings
n <- c(200,500,1000, 2500)
sigma <- c(0.1,0.5)
M <- 100
fam <- c("gaussian", "binomial","poisson")
true.theta <- c(1, 1, 1)/sqrt(3)
params <- rbind(expand_grid(n,sigma, fam = fam[1]),
                expand_grid(n,sigma = 1, fam = fam[-1]))
set.seed(1234)
seeds <- sample.int(5000, M)

### loop on all parameters values
tmp0 <- lapply(1:nrow(params), function(ii, cc){
  print(paste0("----", ii,"/", nrow(params), "----"))
  results <- mclapply(seeds,
                      FUN = function(ss){
                        tryCatch({fit_models(seed = ss, params = params[ii,], cc = cc)},
                                 error = function(err) {return(as.character(err))})
                      }, mc.cores = 4)
  results <- results[which(sapply(results, length)>1)]
  res <- data.frame(c(sapply(results, function(xx) xx$gamnl$single_test), 
                     sapply(results, function(xx) xx$gplsim$single_test),
                     sapply(results, function(xx) xx$gam$single_test)))
  colnames(res) <- "test"
  res[,"Model"] <- factor(rep(c("gamnl", "gplsim", "gam"), each = length(results)))
  res[,"n"] <- factor(paste0("n = ", params[ii, "n"]))
  res[,"cc"] <- factor(cc)
  res[,"sigma"] <- factor(paste0("sigma = ", params[ii, "sigma"]))
  res[,"family"] <- factor(as.vector(params[ii, "fam"]))
  
  return(list(results = res))
}, cc=0) ## end lapply loop

save(tmp0, file = "./results/gplsims_test.rda")

## hist test
test <- do.call("rbind", lapply(tmp0, function(xx) xx$results))

newDF <- as.data.frame(test) %>% 
  pivot_longer(col= !contains(c("Model", "n", "sigma", "fam", "cc")), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))
newDF$Model <- stringr::str_to_title(newDF$Model)
newDF$family <- stringr::str_to_title(newDF$family)
cols <- viridis::viridis(n = 3)

th <- theme(axis.text = element_text(size=15),
            legend.title = element_text(size=17),
            axis.title=element_blank(),
            legend.text = element_text(size=15),
            title = element_text(size=17))
newDF %>% ggplot(aes( sample = column_values, color = Model)) + 
  stat_qq(distribution = stats::qunif) +
  scale_color_manual(values = cols)+
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  facet_grid(cols = vars(n), rows = vars(sigma, family), scales = "free_y")+th+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        strip.text = element_text(size = 17))+xlab("") +
  ggtitle("P-value distribution under H0")
