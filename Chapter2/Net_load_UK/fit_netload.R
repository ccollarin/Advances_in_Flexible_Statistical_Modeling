library(gamFactory)
library(mgcViz)
library(scoringRules)
load("./data/electricity_data.rdata")

df_rows <- dat$Date<"2018-01-01"
df <- dat[df_rows,]
test <- dat[!df_rows,] # test set on 2018 data

## Define formulas
formulas <- list(snest = list(node_nt ~
                                t + dow_RpH_A + School_Hol_A  + 
                                s(doy, bs = "ad", k = 30) + 
                                s(SSRD_t) + s(TP_t)+s(x2ts_t) +
                                s_nest(Wind100_2_Cap, trans = trans_linear(pord = 1)) +
                                s_nest(tempMat, trans = trans_nexpsm())+
                                s_nest(lagMat, trans = trans_linear(pord=1)), 
                              ~ s(doy) + School_Hol_A + storm),
                 naive = list(node_nt ~
                                t + dow_RpH_A + School_Hol_A  + 
                                s(doy, bs = "ad", k = 30) +
                                s(SSRD_t) + s(x2ts_t) + s(Wind100_2_Cap_t) +
                                s(TP_t) + s(node_nt_L1),
                              ~ s(doy) + School_Hol_A  + storm)
)

if (file.exists("fitted_UK.rda")) {
  load("fitted_UK.rda")
  print(fitting_time)
} else {
  fits <- parallel::mclapply(formulas, function(ff){
    start <- format(Sys.time())
    fit <- gam_nl(ff,family = fam_gaussian(),
                  data = df, optimizer = "bfgs")
    time <- difftime(format(Sys.time()),start,  units = "secs")
    return(list(fit = fit, time = time))
  },mc.cores = 2)
  fitting_time <- lapply(fits, function(ff) ff$time)
  fits <- lapply(fits, function(ff) ff$fit)
  names(fits) <- names(formulas)
  print(fitting_time)
  save(fits, fitting_time, file = "fitted_UK.rda")
}

## Summary
summary(fits[[1]])
summary(fits[[2]])
  
## Checks
fits[[1]] <- getViz(fits[[1]], nsim = 100)
fits[[2]] <- getViz(fits[[2]], nsim = 100)

check(fits[[1]])
check(fits[[2]])

nsim <- 20
data.frame(
  logloss = -sapply(fits, function(xx) logLik(xx)),
  crps = sapply(fits, function(ff){
    sum(crps_sample(y = ff$y, 
                    dat = t(sapply(1:nrow(ff$fitted.values), 
                                   function(ii) 
                                     rnorm(n=nsim, mean = ff$fitted.values[ii,1],
                                           sd = sqrt(ff$fitted.values[ii,2]))))))
  }),
  AIC = sapply(fits, AIC),
  test_RMSE = sapply(fits, function(ff) sd(predict(ff, newdata = test)[,1] - test$node_nt)),
  time = unlist(fitting_time)
  )

## Plots
plogis(fits[[1]]$coefficients["s(tempMat).2"])

dd <- data.frame("GSP Group ID" = c("A","B","C","D","E","F","G","P","N","J","H","K","L","M"),
                 Area = c("East England","East Midlands","London","North Wales, Merseyside and Cheshire",
                          "West Midlands","North East England",
                          "North West England","North Scotland","South and Central Scotland",
                          "South East England","Southern England","South Wales National",
                          "South West England","Yorkshire"))
dd[order(dd$GSP.Group.ID),]

### Exponential smooth of temperature ####
th <- theme(axis.text = element_text(size=22),
            axis.title=element_text(size=25), 
            legend.title = element_text(size=22), legend.text = element_text(size=18))
p <- sapply(1:length(fits[[1]]$smooth), function(ii) print(plot(fits[[1]], select = ii)+th))

plot(fits[[1]], select = 6)+ l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciPoly() + l_fitLine() + th + 
  ylab(expression(s(tilde(s)[exp](temp)))) +
  xlab(expression(tilde(s)[exp](temp)))

ss <- which(sapply(fits[[1]]$smooth,"[[", "label") == "s(tempMat)")
si <- fits[[1]]$smooth[[ss]]$xt$si
prange <- (fits[[1]]$smooth[[ss]]$first.para:fits[[1]]$smooth[[ss]]$last.para)[2]
Va <- fits[[1]]$Vp[prange, prange, drop = FALSE]

XX <- test$tempMat
nms <- colnames(XX)

stemp <- expsmooth(y = as.vector( XX[ , which(nms == "y")] ), 
                   Xi =  matrix(as.vector(XX[ , which(nms == "x")] ), ncol=1),
                   beta = si$alpha[2], deriv = 1)

tmp <- data.frame(d0 = stemp$d0,
                  se = sqrt(pmax(0, rowSums((stemp$d1 %*% Va) * stemp$d1))),
                  index = 1:nrow(XX),
                  temp = as.vector(XX[,which(colnames(XX)=="y")]))
tmp$u <- tmp$d0 + qnorm(0.975) * tmp$se
tmp$l <- tmp$d0 - qnorm(0.975) * tmp$se

ggplot(tmp, aes(x = index, y=d0))+
  geom_ribbon(aes(ymin = l, ymax = u), fill = "grey", alpha=0.8) +
  geom_line() +
  geom_line(data = tmp, mapping = aes(x = index, y=temp), linetype = "dashed", color = "darkred")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top") + 
  xlab("Day of the Year")+ylab("Temperature")+th

### wind single index ####
plot(fits[[1]], select = 5) +
  l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciPoly() + l_fitLine()+th + 
  ylab(expression(paste(s(wcap%.%wsp^T, alpha[wsp])))) +
  xlab(expression(paste(wcap%.%wsp^T, alpha[wsp])))

plot(fits[[1]], inner = TRUE, select = 5) + 
  geom_abline(slope=0,intercept = 0, alpha = 0.5) +
  scale_x_discrete(labels = regions)+ th+
  xlab("GSP group") + ylab(expression(alpha[wsp]))


### 
plot(fits[[1]], select = 7) +
  l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciPoly() + l_fitLine()+th + 
  ylab(expression(paste(s(wcap%.%wsp^T, alpha[wsp])))) +
  xlab(expression(paste(wcap%.%wsp^T, alpha[wsp])))

plot(fits[[1]], inner = TRUE, select = 7)+ 
  geom_abline(slope=0,intercept = 0, alpha = 0.5) + 
  scale_x_discrete(labels=gsub("node_nt_L1_","", colnames(dat$lagMat))) + th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

