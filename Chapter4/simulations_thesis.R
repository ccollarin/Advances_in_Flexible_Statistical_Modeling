rm(list = ls())
# install.packages("Boosting4Gam_0.1.0.tar.gz", repos = NULL, type="source")
library(Boosting4Gam)
library(mboost)
source("DeselectBoost.R")
n <- 2000
nsim <- 100
p <- 12
f0 <- function(x,y) (x*y)^2
f1 <- function(x,y) cos(x+y)
f2 <- function(x,y) exp(x-y)
f3 <- function(x) x*(1+x)^3
f4 <- function(x) exp(x/4+1)^2
f5 <- function(x) -x

set.seed(1234)
seeds <- sample.int(5000, nsim)

fun.to.sim <- function(ss, sigma = 2){
  set.seed(ss)
  X <- matrix(rnorm(n*p), ncol=p)
  colnames(X) <- paste0("x", 0:(p-1))
  f <- scale(1/3*f0(X[, "x0"],X[, "x1"]) + 2*f1(X[, "x2"],X[, "x3"]) +
       1/5*f2(X[, "x4"],X[, "x5"]) + 1/16*f3(X[, "x6"]) +
        1/4*f4(X[, "x7"]) + f5(X[, "x8"]), scale = FALSE)

  dat <- as.data.frame(cbind(f+rnorm(n, sd = sigma),X))
  colnames(dat)[1] <- "y"
  val <- dat[(n/2+1):n,]
  dat <- dat[1:(n/2),]
  form <- as.formula(paste("y~",
                           paste0(paste0("s(", colnames(X),")"),collapse = "+")))

  fit <- interactions.sel(form, dat, val , nu.type = 0.01)
  fit.ad <- interactions.sel(form, dat, val, type="l0", nu.type = 0.01)

  out <- list(sigma = sigma)
  out$fixed <- list(step1 = unname(names(fit$fit1$gained.loss)[fit$fit1$gained.loss>0]),
                    step2 = fit$cov1$effects,
                    step3 = unname(names(fit$fit2$gained.loss)[fit$fit2$gained.loss>0]),
                    step4 = fit$cov2$effects,
                    m1 = fit$fit1$m, m2 = fit$fit2$m)
  out$adaptive <- list(step1 = unname(names(fit.ad$fit1$gained.loss)[fit.ad$fit1$gained.loss>0]),
                       step2 = fit.ad$cov1$effects,
                       step3 = unname(names(fit.ad$fit2$gained.loss)[fit.ad$fit2$gained.loss>0]),
                       step4 = fit.ad$cov2$effects,
                       m1 = fit$fit1$m, m2 = fit$fit2$m)

  ### mboost
  out$mboost <- list()
  form.mboost <- as.formula(paste("y~",
                                  paste0(paste0("bbs(", colnames(X),")"),collapse = "+")))

  fit1 <- gamboost(form.mboost, data = dat, boost_control(nu = 0.01, mstop = 1000))
  out$mboost$step1 <- names(coef(fit1))
  out$mboost$step1 <- gsub("bbs", "s", out$mboost$step1)
  gb.des <- DeselectBoost(fit1, fam = Gaussian(), data = dat)
  out$mboost$step2 <- names(coef(gb.des))
  out$mboost$step2 <- gsub("bbs", "s", out$mboost$step2)
  couples <- combn(unique(gsub(".*\\((.*)\\).*", "\\1", names(coef(gb.des)))),2)
  interactions <- paste0(sapply(1:ncol(couples),
                                function(jj) paste0("bspatial(", couples[1,jj], ", ",
                                                    couples[2,jj], ")")), collapse = "+")
  form.mboost2 <-as.formula(paste(all.vars(form.mboost)[1], "~ ",
                                  paste0(c(names(coef(gb.des)),
                                           interactions),  collapse = "+")))
  fit2 <- gamboost(form.mboost2, data = dat, boost_control(nu = 0.01, mstop = 1000))
  gb.des2 <- DeselectBoost(fit2, fam = Gaussian(), data = dat)
  out$mboost$step4 <- names(coef(gb.des2))
  out$mboost$step4 <- gsub("bbs", "s", out$mboost$step4)
  out$mboost$step4 <- gsub("t\\)", "t, k=5\\)", out$mboost$step4)
  out$mboost$step4 <- gsub("bspatial", "ti", out$mboost$step4)

  return(out)
}

library(parallel)
sims <- mclapply(seeds, fun.to.sim, mc.cores = 4)
sims3 <- mclapply(seeds, fun.to.sim, mc.cores = 4, sigma = 3)
sims4 <- mclapply(seeds, fun.to.sim, mc.cores = 4, sigma = 4)
sims5 <- mclapply(seeds, fun.to.sim, mc.cores = 4, sigma = 5)

sims <- sims[sapply(sims, function(ll) length(ll)>1)]
sims3 <- sims3[sapply(sims3, function(ll) length(ll)>1)]
sims4 <- sims4[sapply(sims4, function(ll) length(ll)>1)]
sims5 <- sims5[sapply(sims5, function(ll) length(ll)>1)]
## fix
fixed <- lapply(sims, "[[", "fixed")
fixed3 <- lapply(sims3, "[[", "fixed")
fixed4 <- lapply(sims4, "[[", "fixed")
fixed5 <- lapply(sims5, "[[", "fixed")
post.proc <- function(lsims, step, sigma, truth){
  effs <- unique(unlist(lapply(lsims, "[[", step)))
  if(step<3) truth <- truth[grep("s\\(|bbs\\(", truth)]

  nwrong <- sapply(lsims, function(ll) sum(!gsub("bbs\\(", "s\\(", ll[[step]])%in%truth))
  nwright <- sapply(lsims, function(ll) sum(gsub("bbs\\(", "s\\(", ll[[step]])%in%truth))
  meff <- rep(0, length = length(effs))
  names(meff) <- effs
  for(ss in lsims){
    meff[ss[[step]]] <- meff[ss[[step]]] +1
  }
  return(list(mat = c(sigma = sigma, meff),
              nwrong = cbind(sigma = sigma, n = nwrong),
              nwright = cbind(sigma = sigma, n = nwright)))
}

true1 <- paste0("s(x", c(0:8), ")")
true4 <- c(true1, "ti(x0,x1)", "ti(x2,x3)","ti(x4,x5)")
step.fix <- list(lapply(paste0("step", 1:4), function(st) post.proc(fixed, st,
                                                                    sigma = 2, truth = true4)),
                 lapply(paste0("step", 1:4), function(st) post.proc(fixed3, st,
                                                                    sigma = 3, truth = true4)),
                 lapply(paste0("step", 1:4), function(st) post.proc(fixed4, st,
                                                                    sigma = 4, truth = true4)),
                 lapply(paste0("step", 1:4), function(st) post.proc(fixed5, st,
                                                                    sigma = 5, truth = true4)))

mat.fix <- lapply(1:length(step.fix), function(ii){
  ee <- lapply(lapply(step.fix, "[[", ii), "[[", "mat")
  mat <- matrix(0, ncol = length(unique(names(unlist(ee)))), nrow=length(step.fix))
  colnames(mat) <- unique(names(unlist(ee)))
  for(ll in 1:length(step.fix)){
    mat[ll, names(ee[[ll]])] <- ee[[ll]]
  }
  return(mat)
} )

library(ggplot2)

th <- theme(axis.text = element_text(size=18),
            axis.title=element_text(size=18),
            legend.title = element_text(size=18),
            legend.text = element_text(size=18),
            title = element_text(size=20),
            strip.text = element_text(size = 18))
# Creazione della heatmap
mat <- (mat.fix[[2]][,-1])[,paste0("s(x", c(0:8, 10), ")")]
cols <- rep("black", ncol(mat))
cols[colnames(mat)%in%true1] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat)
p1.fix <- ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +  theme_minimal() +  labs(x = "", y = "", fill = "") +
  geom_text(aes(label = round(value, 2)), color = "white")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = cols))+
  ylab("LASSO")+th

mat <- (mat.fix[[4]][,-1])[,sort(colnames(mat.fix[[4]])[-1])]
cols <- rep("black", ncol(mat))
cols[colnames(mat)%in%true4] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat)
# Creazione della heatmap
p4.fix <- ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "", y = "", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = cols))+
  ylab("LASSO")+th

################
.mboost <- lapply(sims, "[[", "mboost")
.mboost3 <- lapply(sims3, "[[", "mboost")
.mboost4 <- lapply(sims4, "[[", "mboost")
.mboost5 <- lapply(sims5, "[[", "mboost")
step.mb <- list(lapply(paste0("step", 1:4),
                       function(st) post.proc(.mboost, st, sigma = 2, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(.mboost3, st, sigma = 3, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(.mboost4, st, sigma = 4, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(.mboost5, st, sigma = 5, truth = true4)))
mat.mb <- lapply(1:length(step.mb), function(ii){
  ee <- lapply(lapply(step.mb, "[[", ii), "[[", "mat")
  mat <- matrix(0, ncol = length(unique(names(unlist(ee)))), nrow=length(step.mb))
  colnames(mat) <- unique(names(unlist(ee)))
  for(ll in 1:length(step.mb)){
    mat[ll, names(ee[[ll]])] <- ee[[ll]]
  }
  return(mat)
} )

w.mb <- lapply(1:length(step.mb), function(ii){
  ee <- cbind(do.call("rbind",
                      lapply(lapply(step.mb, "[[", ii), "[[", "nwrong")),
              wright = do.call("rbind",
                               lapply(lapply(step.mb, "[[", ii), "[[", "nwright"))[,-1])
  cbind(model = "mboost",as.data.frame(ee))})


# Creazione della heatmap
mat <- mat.mb[[2]][,-1]
mat <- mat[,sort(colnames(mat))]
cols <- rep("black", ncol(mat))
cols[colnames(mat)%in%true1] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat[,sort(colnames(mat))])
mat.long$Var1 <- mat.long$Var1+1
p1.mb <- ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +  theme_minimal() +  labs(x = "", y = "", fill = "") +
  geom_text(aes(label = round(value, 2)), color = "white")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = cols))+
  ylab("mboost")+th

mat <- (mat.mb[[4]][,-1])[,sort(colnames(mat.mb[[4]])[-1])]
cols <- rep("black", ncol(mat))
true4.mb <- c(true1, "ti(x0, x1)", "ti(x2, x3)","ti(x4, x5)")
cols[colnames(mat)%in%true4.mb] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat)
mat.long$Var1 <- mat.long$Var1+1
# Creazione della heatmap
p4.mb<-ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + theme_minimal() +labs(x = "", y = "", fill = "") +th+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.1, hjust = 1, colour = cols))+
  ylab("mboost")

################
adapt <- lapply(sims, "[[", "adaptive")
adapt3 <- lapply(sims3, "[[", "adaptive")
adapt4 <- lapply(sims4, "[[", "adaptive")
adapt5 <- lapply(sims5, "[[", "adaptive")
step.ada <- list(lapply(paste0("step", 1:4),
                        function(st) post.proc(adapt, st, sigma = 2, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(adapt3, st, sigma = 3, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(adapt4, st, sigma = 4, truth = true4)),
                 lapply(paste0("step", 1:4),
                        function(st) post.proc(adapt5, st, sigma = 5, truth = true4)))
mat.ada <- lapply(1:length(step.ada), function(ii){
  ee <- lapply(lapply(step.ada, "[[", ii), "[[", "mat")
  mat <- matrix(0, ncol = length(unique(names(unlist(ee)))), nrow=length(step.ada))
  colnames(mat) <- unique(names(unlist(ee)))
  for(ll in 1:length(step.ada)){
    mat[ll, names(ee[[ll]])] <- ee[[ll]]
  }
  return(mat)
} )

# Creazione della heatmap
mat <- mat.ada[[2]][,-1]
mat <- mat[,sort(colnames(mat))]
cols <- rep("black", ncol(mat))
cols[colnames(mat)%in%true1] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat)
mat.long$Var1 <- mat.long$Var1+1
p1.ada <- ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +  theme_minimal() +  labs(x = "", y = "", fill = "") +
  geom_text(aes(label = round(value, 2)), color = "white")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = cols))+
  ylab("L0")+th

mat <- (mat.ada[[4]][,-1])[,sort(colnames(mat.ada[[4]])[-1])]
cols <- rep("black", ncol(mat))
cols[colnames(mat)%in%true4] <- "blue"
names(cols) <- colnames(mat)
mat.long <- reshape2::melt(mat)
mat.long$Var1 <- mat.long$Var1+1
# Creazione della heatmap
p4.ada <- ggplot(data = mat.long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + theme_minimal() +labs(x = "", y = "", fill = "") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.1, hjust = 1, colour = cols))+
  ylab("L0")+th

w.l0 <- lapply(1:length(step.ada), function(ii){
  ee <- cbind(do.call("rbind",
                lapply(lapply(step.ada, "[[", ii), "[[", "nwrong")),
        wright = do.call("rbind",
          lapply(lapply(step.ada, "[[", ii), "[[", "nwright"))[,-1])
  cbind(model = "L0",as.data.frame(ee))})


p1 <- gridExtra::grid.arrange(p1.ada, p1.mb)
p4 <- gridExtra::grid.arrange(p4.ada, p4.mb)

b1 <- ggplot(rbind(w.mb[[2]], w.l0[[2]]), aes(x=factor(sigma), y=n))+
  geom_boxplot()+facet_wrap(~factor(model)) + xlab(expression(sigma))+th+
  ggtitle("Step 2")+ylab("n false positives")
b2 <- ggplot(rbind(w.mb[[4]], w.l0[[4]]), aes(x=factor(sigma), y=n))+
  geom_boxplot()+facet_wrap(~factor(model))+th+ggtitle("Step 4")+
  xlab(expression(sigma))+ylab("n false positives")
pbox <- gridExtra::grid.arrange(b1,b2)

b1 <- ggplot(rbind(w.mb[[2]], w.l0[[2]]), aes(x=factor(sigma), y=9 - wright))+
  geom_boxplot()+facet_wrap(~factor(model)) + xlab(expression(sigma))+th+
  ggtitle("Step 2")+ylab("n false negatives")
b2 <- ggplot(rbind(w.mb[[4]], w.l0[[4]]), aes(x=factor(sigma), y=12 - wright))+
  geom_boxplot()+facet_wrap(~factor(model))+th+ggtitle("Step 4")+
  xlab(expression(sigma))+ylab("n false negatives")+ylim(0,13)
(pbox <- gridExtra::grid.arrange(b1,b2))
