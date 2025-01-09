#####################################
# code to perform simulations with gamFactory
#####################################
#devtools::install_github("mfasiolo/mgcViz")
#devtools::install_github("mfasiolo/gamfactory", ref = "dev")

library(gamFactory)
set.seed(1234)
n <- 1000
cc <- 1
nsim <- 10
ncores <- 4
seeds <- sample.int(nsim*10, size=nsim)

sim_gamnl <- function(ss, seeds, cc = c(1,1,1)){
  set.seed(seeds[ss])
  print(ss)
  y <- list()
  ### standard component
  dat <- data.frame(x1 = rnorm(n))
  y[["std"]] <- scale(-0.5*dat$x1^4+2*dat$x1^2)
  
  ### single index 
  p <- 8
  dat$Xsi <- matrix(runif(p * n, -4, 4), n, p)
  aa.si <- c(2.11,  1.39, -1.22, rep(0,p-3))
  aa.si <- aa.si/sum(aa.si^2)
  z.si <- dat$Xsi %*% aa.si
  y[["si"]] <- scale(exp(z.si)+2*sin(z.si*2))*cc[1]
  
  ### exponential smooth
  tim <- seq(0, 2*pi, length.out = n)
  temp <- sin(tim) + rnorm(n, 0, 0.5)
  Xi <- cbind(1, tim, tim^2)
  Xi[,-1] <- scale(Xi[,-1], scale=FALSE)
  aa.expsm <- c(2, 0.1, 0.1)
  z.exp <- expsmooth(y = temp, Xi = Xi, beta = aa.expsm)$d0
  dat$Xexp <- cbind(temp, Xi)
  colnames(dat[,"Xexp"]) <- c("y", rep("x", ncol(Xi)))
  y[["expsm"]] <- scale(z.exp + sin(3 * z.exp))*cc[1]
  
  ### kernel smooth
  n0 <- 50
  X0 <- cbind(runif(n0, -1, 1), runif(n0, -4, 4))
  X <- cbind(runif(n, -1, 1), runif(n, -4, 4))
  dist <- lapply(1:ncol(X), function(dd){
    t(sapply(1:nrow(X), function(ii) (X[ii, dd] - X0[ , dd])^2 ))
  })
  dist <- dist[[1]] + dist[[2]]
  trueF <- function(x) 3 * x[ , 1] + x[ , 2]^2
  temp <- t(sapply(1:n, function(ii) trueF(X0)))
  
  dat$Xks <- cbind(temp,dist)
  colnames(dat$Xks) <- c(rep("y", nrow(X0)), rep("d1", nrow(X0)))
  z.mgks <- trueF(X)
  y[["mgks"]] <- scale(z.mgks*cos(z.mgks/2))*cc[1]
  
  dat$y <- Reduce("+", y) + rnorm(n)
  
  fit <- tryCatch(gam_nl(list(y~s(x1) +
                  s_nest(Xsi, k = 16, trans = trans_linear())+
                  s_nest(Xexp, k = 16, trans = trans_nexpsm())+
                  s_nest(Xks, k = 16, trans = trans_mgks()),
                  ~1), data = dat, family = fam_gaussian(),
                control = list(trace = FALSE), optimizer = "bfgs"),
  warning = function(war) {return(as.character(war))},
  error = function(err) {return(as.character(err))})

  if(length(fit)==1){
    return(fit)
  } else{
    out <- list()
    out$terms <- predict(fit, type = "terms")
    out$summary <- summary(fit)$s.table[,"p-value"]
    out$si <- list(alpha = drop(fit$smooth[[2]]$xt$si$B %*% (fit$smooth[[2]]$xt$si$alpha + fit$smooth[[2]]$xt$si$a0)))
    out$si$xsi <- fit$smooth[[2]]$xt$si$X %*% out$si$alpha
    
    out$smexp <- list(hatw = plogis(fit$smooth[[3]]$xt$si$X %*% 
                      fit$smooth[[3]]$xt$si$alpha[-1]),
                      truew = plogis(dat$Xexp[,-1] %*% aa.expsm))
    out$smexp$xexp <- exp(fit$smooth[[3]]$xt$si$alpha[1]) *
      (expsmooth(y = fit$smooth[[3]]$xt$si$x,
        Xi = fit$smooth[[3]]$xt$si$X,
        beta = fit$smooth[[3]]$xt$si$alpha[-1])$d0-fit$smooth[[3]]$xt$si$xm)
    
    out$smks <- list(xmgks = drop(exp(fit$smooth[[4]]$xt$si$alpha[1])*(
      mgks(y = fit$smooth[[4]]$xt$si$x,
                 dist = fit$smooth[[4]]$xt$si$dist,
                 beta = fit$smooth[[4]]$xt$si$alpha[-1], deriv=0)$d0-
        fit$smooth[[4]]$xt$si$xm)),
      true = z.mgks)
    if(ss==1){
      out$truth <- list(y = y, 
                        x = list(si = scale(z.si), 
                                  expsm= scale(z.exp), 
                                  mgks = scale(z.mgks)))
    }
    return(out)
  }
}

sims <- parallel::mclapply(1:nsim, sim_gamnl, seeds = seeds)

# save(sims, file="./results/full_sims.rda")

sims <- parallel::mclapply(1:nsim, sim_gamnl, seeds = seeds,cc=rep(0,3))

# save(sims, file="./results/full_sims_test0.rda")


## hist test
library(tidyverse)
if(any(sapply(sims, length)==1)) table(sims[sapply(sims, length)==1])
sims <- sims[sapply(sims, length)!=1]
test <- do.call("rbind", lapply(sims, function(xx) xx$summary))

newDF <- as.data.frame(test) %>% 
  pivot_longer(col= contains("s_nest"), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))

cols <- viridis::viridis(n = 3)

newDF %>% ggplot(aes(x = column_values, fill=column_names)) +
  geom_histogram(position='dodge') +
  scale_fill_manual(values = cols)+
  facet_grid(cols = vars(column_names))+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        axis.title = element_blank())+xlab("") +
  ggtitle("P-value distribution under H0")

newDF %>% ggplot(aes( sample = column_values)) + 
  stat_qq(distribution = stats::qunif) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)+
  facet_grid(cols = vars(column_names))+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        axis.title = element_blank())+xlab("") +
  ggtitle("P-value distribution under H0")

### alpha boxplot
alpha <- data.frame(do.call("rbind", lapply(sims, function(xx) xx$si$alpha/sqrt(sum(xx$si$alpha^2)))))
aa.si <- c(2.11,  1.39, -1.22, rep(0,ncol(alpha)-3))
aa.si <- aa.si/sqrt(sum(aa.si^2))
newDF <- as.data.frame(alpha) %>% 
  pivot_longer(col= contains("X"), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))

cols <- viridis::viridis(n = 3)

(p1 <- newDF %>% ggplot(aes(x = column_names, y = column_values)) +
    geom_boxplot() + geom_abline(intercept = 0, slope=0, color = "darkgrey") +
    geom_point(data.frame(y = aa.si,x = factor(paste0("X", 1:length(aa.si)))),
               mapping = aes(x=x,y=y), color="red")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(), legend.position = "none")+xlab("") +
    ggtitle("Estimated single index parameters"))

effects <- do.call("rbind", lapply(1:length(sims), 
                                   function(ii) data.frame(x=sims[[ii]]$si$xsi, 
                                                           y = sims[[ii]]$terms[,"s(Xsi)"],
                                                           sim = ii)))
(p1 <- effects %>% 
    ggplot(aes(y = y, x = x, fill = factor(sim))) +
    geom_line(color="grey") +
    geom_line(data.frame(x=sims[[1]]$truth$x$si, 
                         y=sims[[1]]$truth$y$si -
                           sims[[1]]$truth$y$si[which.min(abs(sims[[1]]$truth$x$si))],
                         sim = 1), 
              mapping = aes(x=x,y=y, fill = ), color= "red")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(), legend.position = "none")+xlab("") +
    ggtitle("Single index effect"))
### nexpsm 
effects <- do.call("rbind", lapply(1:length(sims), 
                                   function(ii) data.frame(x=sims[[ii]]$smexp$xexp, 
                                                           y = sims[[ii]]$terms[,"s(Xexp)"],
                                                           sim = ii)))
(p2 <- effects %>% 
    ggplot(aes(y = y, x = x, fill = factor(sim))) +
    geom_line(color="grey") + 
    geom_line(data.frame(x=sims[[1]]$truth$x$expsm, 
                         y=sims[[1]]$truth$y$expsm -
                           sims[[1]]$truth$y$expsm[which.min(abs(sims[[1]]$truth$x$expsm))],
                         sim = 1), 
              mapping = aes(x=x,y=y, fill = ), color= "red")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(), legend.position = "none")+xlab("") +
    ggtitle("Exponential smoothing effect"))
### mgks 
effects <- do.call("rbind", lapply(1:length(sims), 
                                   function(ii) data.frame(x=sims[[ii]]$smks$xmgks, 
                                                           y = sims[[ii]]$terms[,"s(Xks)"],
                                                           sim = ii)))
(p3 <- effects %>% 
    ggplot(aes(y = y, x = x, fill = factor(sim))) +
    geom_line(color="grey") + 
    geom_line(data.frame(x=sims[[1]]$truth$x$mgks, 
                         y=sims[[1]]$truth$y$mgks -
                           sims[[1]]$truth$y$mgks[which.min(abs(sims[[1]]$truth$x$mgks))],
                         sim = 1), 
              mapping = aes(x=x,y=y, fill = ), color= "red")+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(), legend.position = "none")+xlab("") +
    ggtitle("Kernel smoothing effect"))


ggsave(gridExtra::grid.arrange(p1,p2,p3, ncol=3),
        filename = "./figures/sim_full_effects.pdf",
       width = 13, height = 5, units = "in")
