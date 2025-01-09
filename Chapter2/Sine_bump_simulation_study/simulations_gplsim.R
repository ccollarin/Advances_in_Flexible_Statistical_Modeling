########################################
#  Simulations to compare results from gplsim package and gamFactory
########################################
rm(list=ls())
library(gplsim)
library(cgaim)
library(gamFactory)
library(parallel)
library(ggplot2)
library(tidyverse)
source("./simulations_gplsim_functions.R")

# Gaussian family
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
tmp <- lapply(1:nrow(params), function(ii){
  print(paste0("----", ii,"/", nrow(params), "----"))
  results <- mclapply(seeds, fit_models, params = params[ii,], mc.cores = 4)
  nmod <- names(results[[1]])[names(results[[1]])%in%c("gamnl", "gplsim", "cgaim")]
  res <- data.frame(do.call("rbind", 
            lapply(nmod, function(nn) t(sapply(results, function(xx) xx[[nn]]$coefs)))))
  colnames(res) <- c("theta1", "theta2", "theta3", "gamma")
  res[,"Model"] <- factor(rep(nmod, each = M))
  res[,"n"] <- factor(paste0("n = ", params[ii, "n"]))
  res[,"sigma"] <- factor(paste0("sigma = ", params[ii, "sigma"]))
  res[,"family"] <- factor(as.vector(params[ii, "fam"]))
  res[,"mse"] <- do.call("c", 
                         lapply(nmod, function(nn) t(sapply(results, function(xx) xx[[nn]]$mse))))
  res[,"time"] <- do.call("c", 
                          lapply(nmod, function(nn) t(sapply(results, function(xx) xx[[nn]]$time))))

  coefs2 <- sapply(results, function(xx) c(xx$gamnl$coefs2))
  return(list(results = res,
              coefs2 = coefs2))
}) ## end lapply loop

## boxplot coefficients
coefs1 <- do.call("rbind", lapply(tmp, function(xx) xx$results[, c("theta1", "theta2", "theta3", "gamma", 
                                                                   "Model", "n", "sigma", "family")]))
coefs1$theta1 <- coefs1$theta1 - true.theta[1]
coefs1$theta2 <- coefs1$theta2 - true.theta[2]
coefs1$theta3 <- coefs1$theta3 - true.theta[3]
coefs1$gamma <- coefs1$gamma - 0.3
coefs1$family <- stringr::str_to_title(coefs1$family)

newDF <- as.data.frame(coefs1) %>% 
  pivot_longer(col= !contains(c("Model", "n", "sigma", "fam")), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))


cols <- viridis::viridis(n = 3)
th <- theme(axis.text = element_text(size=15),
      axis.title=element_blank(),
      legend.title = element_text(size=17),
      legend.text = element_text(size=15),
      title = element_text(size=17),
      strip.text = element_text(size = 17))

(p1 <- newDF %>% ggplot(aes(x = n, y = column_values, fill = Model)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0, color = "darkgrey") +
    scale_fill_manual(values = cols)+
    facet_grid(cols = vars(column_names), rows = vars(sigma, family), scales = "free_y")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(),)+th)

## boxplot mse
mse1 <- do.call("rbind", lapply(tmp, function(xx) xx$results[,c("mse", "Model", "n", "sigma", "family")]))
newDF <- as.data.frame(mse1) %>% 
  pivot_longer(col= !contains(c("Model","sigma","n", "family")), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))

(p2 <- newDF %>% ggplot(aes(x = n, y = column_values, fill = Model)) +
    geom_boxplot() +
    scale_fill_manual(values = cols)+
    # facet_wrap(facets = vars(sigma, family), scales = "free_y")+
    facet_grid(cols = vars(column_names), rows = vars(family, sigma), scales = "free_y")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank())+xlab("")+th
)

## boxplot time
time1 <- do.call("rbind", lapply(tmp, function(xx) xx$res[c("time", "Model", "n", "sigma", "family")]))
newDF <- as.data.frame(time1) %>% 
  pivot_longer(col= !contains(c("Model","sigma","n", "family")), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))
newDF$family <- stringr::str_to_title(newDF$family)

scaleFUN <- function(x) sprintf("%.2f", x)
(p4 <- newDF %>% ggplot(aes(x = n, y = column_values, fill = Model)) +
    geom_boxplot() +
    scale_fill_manual(values = cols)+
    scale_y_continuous(trans='log',labels = scaleFUN )+
    facet_grid(cols = vars(column_names), rows = vars(sigma, family), scales = "free_y")+
    theme(axis.text.x = element_text(angle = 60, hjust=1))+th+
    xlab("")+ylab("Seconds")
)

###################################
# summarizing tables
###################################
library(dplyr)
coefs1$theta1 <- coefs1$theta1 + true.theta[1]
coefs1$theta2 <- coefs1$theta2 + true.theta[2]
coefs1$theta3 <- coefs1$theta3 + true.theta[3]
coefs1$gamma <- coefs1$gamma + 0.3
coefs1 %>% group_by(sigma, n, Model, family) %>% 
  summarise(across(!contains(c("Model","sigma","n")),
                   list(mean = mean, sd = sd)))
mse1 %>% group_by(sigma, n, Model, family) %>% 
  summarise(across(!contains(c("Model","sigma","n", "family")),
                   list(mean = mean, sd = sd))) %>%
  arrange(family, sigma, n, Model)

time1 %>% group_by(sigma, n, Model, family) %>% 
  summarise(across(!contains(c("Model","sigma","n", "family")),
                   list(mean = mean, sd = sd)))

