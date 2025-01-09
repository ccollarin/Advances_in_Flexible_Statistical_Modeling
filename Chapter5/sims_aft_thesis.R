rm(list = ls())
# install.packages("aftNested_0.1.0.tar.gz", repos = NULL, type="source")
library(aftNested)

N <- 2000
nsim <- 100

fun.sim <- function(ss){
  set.seed(ss)
  out <- list()
  aa <- c(1,-1,1,1)
  dat <- simulate.data(N, smooth.term=2)
  p <- 4

  sample_val <- sample.int(N, N/2)
  validation <- dat[sample_val,]
  dat <- dat[-sample_val,]

  fit2 <- tryCatch({aft_nl(logT ~ s(X5, k=10) + X1 + X2 + X3+X4,
                 outer.args = list(k=15),
                 data = dat, val = validation, unif=FALSE,
                 censored = "Event", trace = TRUE, ps = TRUE)},
                 error = function(err) {return(as.character(err))})

  if(is.character(fit2)){
    out$fit2 <- fit2
    return(out)
  }

  f2 <- get.si(fit2$smooth, fit2$coefficients)

  out$lin.coef <- f2$alpha[1:p+1]
  out$sm <- cbind(y = -fit2$smooth$xt$si$X[,-(1:(p+1))]%*%(f2$alpha[-(1:(p+1))]/fit2$coefficients[1]),
                  x = dat$X5,
                  seed = ss)
  out$s <- cbind(s = 1-exp(f2$s$s),
                 x = f2$xba,
                 seed = ss)
  out$lambda <- fit2$lambda
  out$edf <- c(outer = sum(fit2$edf[-(1:length(f2$alpha))]),
               inner = sum(fit2$edf[fit2$smooth$xt$si$fplp[[1]]]))
  return(out)
}
set.seed(1234)
seeds <- sample.int(5000, nsim)
library(parallel)
sims <- mclapply(seeds, fun.sim, mc.cores = 4)
save(sims, file = "sims_aft.rda")

sims <- sims[sapply(sims, function(ss) length(ss)>1)]

lambda <- t(sapply(sims, "[[", "lambda"))
library(tidyverse)
th <- theme(axis.text = element_text(size=15),
            legend.title = element_text(size=17),
            legend.text = element_text(size=15),
            title = element_text(size=17),
            strip.text = element_text(size = 17))
ggplot(pivot_longer(as.data.frame(lambda), cols=c(outer, inner)),
       aes(x=factor(name), y=exp(value)))+
  geom_boxplot()+
  xlab("")+ylab(expression(lambda))+th

lin.coef <- t(sapply(sims, "[[", "lin.coef"))
colnames(lin.coef) <- paste0("a", 1:4)
library(tidyverse)
ggplot(pivot_longer(as.data.frame(lin.coef),
                    cols=colnames(lin.coef)), aes(x=factor(name), y=value))+
  geom_boxplot()+
  xlab("")+ylab("")+geom_point(data.frame(name=factor(colnames(lin.coef)),
                                          value = c(1,-1,1,1)),
                               mapping=aes(x=name, y=value), color ="red", size = 2)+th

# smooth
sms <- do.call("rbind", lapply(sims, "[[", "sm"))
colnames(sms)[1] <- "y"
dat <- simulate.data(N/2, smooth.term=2)
ggplot(sms, aes(x=x, y=y, color = factor(seed)))+
  geom_line()+
  xlab(expression(X[5]))+ylab("")+
  geom_line(data.frame(x=dat$X5,y=dat$X5^2-mean(dat$X5^2), seed=1),
            mapping=aes(x=x, y=y, color=seed), color ="red", linewidth = 2)+
  theme(legend.position = "none")+
  scale_color_viridis_d()+th

# survival
sms <- do.call("rbind", lapply(sims, "[[", "s"))
ggplot(sms, aes(x=x, y=1-s, color = factor(seed)))+
  geom_line()+
  xlab("Estimated log time")+ylab("Survival function")+
  geom_line(data.frame(x=log(dat$Texp),y=1-rank(dat$Texp)/length(dat$Texp)),
            mapping=aes(x=x, y=y), color ="red", linewidth = 2)+
  theme(legend.position = "none")+
  scale_color_viridis_d()+th

