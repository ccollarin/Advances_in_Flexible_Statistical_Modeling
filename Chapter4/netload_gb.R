rm(list=ls())
# install.packages("Boosting4Gam_0.1.0.tar.gz", repos = NULL, type="source")
library(Boosting4Gam)
library(mboost)
source("DeselectBoost.R")

fun.fit <- function(dat){

  cat("\n#######################################\n")
  dat <- dat[dat$clock_hour==12, ]

  val_set <- dat$Date>as.Date("2016-12-31") & dat$Date<=as.Date("2017-12-31")
  test_set <- dat$Date>as.Date("2017-12-31")
  val <- dat[val_set,]
  test <- dat[test_set, ]
  dat <- dat[!(val_set | test_set),]
  form <- node_n~ s(doy, bs= "cc", k=20) + s(t, k=5)+
    s(node_L1) + s(SSRD) + s(Wind100_2_Cap) + s(x2T) + s(TP)

  out <- list()
  out$fit <- interactions.sel(form, dat, val, nu.type = 0.01)
  out$fit.ad <- interactions.sel(form, dat, val, type="l0", nu.type = 0.1)

  ### mboost
  out$mboost <- list()
  form.mboost <- node_n~ bbs(doy) + bbs(t)+ bbs(node_L1) + bbs(SSRD) + bbs(Wind100_2_Cap) + bbs(x2T) + bbs(TP)

  out$mboost$fit1 <- gamboost(form.mboost, data = dat, boost_control(nu = 0.01, mstop = 1000))
  gb.des <- DeselectBoost(out$mboost$fit1, fam = Gaussian(), data = dat)
  couples <- combn(unique(gsub(".*\\((.*)\\).*", "\\1", names(coef(gb.des)))),2)
  interactions <- paste0(sapply(1:ncol(couples),
                                function(jj) paste0("bspatial(", couples[1,jj], ", ",
                                                    couples[2,jj], ")")), collapse = "+")
  form.mboost2 <-as.formula(paste(all.vars(form.mboost)[1], "~ ",
                                  paste0(c(names(coef(gb.des)),
                                           interactions),  collapse = "+")))
  out$mboost$fit2 <- gamboost(form.mboost2, data = dat, boost_control(nu = 0.01, mstop = 1000))
  gb.des2 <- DeselectBoost(out$mboost$fit2, fam = Gaussian(), data = dat)
  out$mboost$cov2 <- names(coef(gb.des2))
  out$mboost$cov2 <- gsub("bbs", "s", out$mboost$cov2)
  out$mboost$cov2 <- gsub("t\\)", "t, k=5\\)", out$mboost$cov2)
  out$mboost$cov2 <- gsub("bspatial", "te", out$mboost$cov2)
  out$mboost$final.formula <- as.formula(paste0("node_n~",
                                    paste0(out$mboost$cov2, collapse = "+")))
  out$mboost$cov2 <- gsub("te\\(", "ti\\(", out$mboost$cov2)

  out$gam <- gam(form,
                 data = rbind(dat, val))
  out$fit.gam <- gam(out$fit$final.formula,
                 data = dat)
  ff <- as.formula(paste0("node_n~", as.character(out$fit.ad$final.formula)[3]))
  out$fit.ad.gam <- gam(ff,
                    data = dat)
  out$fit.mboost.gam <- gam(out$mboost$final.formula,
                    data = dat)

  out$yg <- predict(out$gam, test)
  out$y1 <- predict(out$fit.gam, test)
  out$y.ad <- predict(out$fit.ad.gam, test)
  out$y.mboost <- predict(out$fit.mboost.gam, test)

  out$loss <- c(sum((out$yg-test$node_n)^2),
                sum((out$y1-test$node_n)^2),
                sum((out$y.ad-test$node_n)^2),
                sum((out$y.mboost-test$node_n)^2))
  return(out)
}

fits <- lapply(netload, fun.fit)

eff.fit <- sapply(fits, function(ff) unname(ff$fit$cov2$effects))
eff.ad <- sapply(fits, function(ff) unname(ff$fit.ad$cov2$effects))
eff.mboost <- sapply(fits, function(ff) unname(ff$mboost$cov2))

eff.sel <- matrix(0, nrow = 8, ncol=length(unique(unlist(c(eff.ad, eff.mboost)))))
colnames(eff.sel) <- sort(unique(unlist(c(eff.ad, eff.mboost))))
nn <- c("East England", "London", "South West England", "North Scotland")
rownames(eff.sel) <- c(paste0(nn, "_0"),
                       paste(nn, "_mboost"))

for(ll in 1:8){ #12
  tmp <- c(eff.ad, eff.mboost)[[ll]]
  eff.sel[ll, tmp] <- 1
}
eff.sel <- eff.sel[sort(rownames(eff.sel)),]

library(ggplot2)
# Creazione della heatmap
mat.long <- reshape2::melt(t(eff.sel))
th <- theme(axis.text = element_text(size=18),
            axis.title=element_text(size=18),
            legend.title = element_text(size=22),
            legend.text = element_text(size=18),
            title = element_text(size=22))
ggplot(data = mat.long, aes(y = Var2, x = Var1, fill = factor(value))) +
  geom_tile() +  theme_minimal() +  labs(x = "", y = "", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_viridis_d(labels=c("non selected", "selected"))+th

xtable::xtable(sapply(fits, "[[", "loss"))
