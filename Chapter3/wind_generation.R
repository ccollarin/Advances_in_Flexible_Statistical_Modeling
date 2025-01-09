rm(list=ls())
library(gamFactory)
library(mgcViz)
load("./wind_gen.rda")

## transform the response and its lags
wind.gen$yt <- log((wind.gen$h2pot)/(max(wind.gen$h2pot)+1e-8 -wind.gen$h2pot)+1e-8)
wind.gen$yt_h24 <- log((wind.gen$lag_h24)/(max(wind.gen$lag_h24)+1e-8 -wind.gen$lag_h24)+1e-8)
wind.gen$yt_h1 <- log((wind.gen$lag_h1)/(max(wind.gen$lag_h1)+1e-8 -wind.gen$lag_h1)+1e-8)

## prepare metrices for nested effects
wind.gen$wind_80 <- as.matrix(cbind(log(wind.gen[, grep("wind_80m_", colnames(wind.gen))]+1), 
                                    wind.gen$dists))
colnames(wind.gen$wind_80) <- c(rep("y", 10), rep("d1", 10))
wind.gen$wdir_80 <- as.matrix(wind.gen[, grep("wdir_80m_", colnames(wind.gen))[-9]])

dat <- wind.gen

dat0 <- dat[dat$h2pot>0, ]
plot(density(dat0$yt))

fit <- gam_nl(list(yt~impianto + month + s(hour) + s(yt_h1) + 
                 s_nest(wind_80, trans = trans_mgks()) +
                 s_nest(wdir_80, trans = trans_linear())
                  , ~1,~1,~1),
               data = dat0, family = fam_shash(), optimizer="bfgs", control=list(trace=TRUE))

fit <- getViz(fit)
summary(fit)

print(plot(fit), pages=1)
print(plot(fit, inner = TRUE, select = 3:4), pages=1)

check(fit)
