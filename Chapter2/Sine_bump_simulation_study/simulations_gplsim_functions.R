### generate data from sine-bump function as described in Carrol et al. (1997)
generate_data.gplsim <- function (n, family = "gaussian", scale,
                                  true.theta = c(1, 1, 1)/sqrt(3), 
                                  c1 = 0.3912, c2 = 1.3409, rho = 0.3, cc=1){
  X <- matrix(stats::runif(length(true.theta) * n), ncol = length(true.theta))
  U <- X %*% true.theta
  fU <- sin((U - c1) * pi/(c2 - c1))*cc
  Z <- 1 - c(1:n)%%2
  q <- as.vector(fU + rho * Z)
  
  if (family == "gaussian") {
    y <- q + rnorm(n, 0, scale)
  }else if(family=="binomial"){
    py <- plogis(q)
    y <- rbinom(length(q), size=1, prob=py)
  }else if(family=="poisson"){
    py <- exp(q)
    y <- rpois(length(q),py)
  }else if(family=="gamma"){
    mu <- q +1
    alpha <- exp(q)
    y <- rgamma(length(q), shape=alpha, scale=mu/alpha)
  }
  
  dat <- data.frame(y = y, ax = U)
  dat$X <- X # single-index predictors
  dat$Z <- Z # partially linear predictors
  
  return(dat)
} ## end generate_data.gplsim

### fit gam_nl and gplsim models 
fit_models <- function(seed, params, cc=1){
  set.seed(seed)
  nn <- as.numeric(params[1])
  ss <- as.numeric(params[2])
  fam <- as.character(params[3])
  # built dataset to fit models
  dat <- generate_data.gplsim(nn, family=fam, scale = ss, cc=cc)
  test <- generate_data.gplsim(nn, family=fam, scale = ss, cc=cc)
  
  form_gamnl <- list(y ~ s_nest(X, k = 10, trans = trans_linear(), m = c(4, 2)) + as.factor(Z))
  if(grepl("gauss", fam)){
    form_gamnl[[2]] <- ~1
    fam_gplsim <- gaussian
    fam_gamnl <- fam_gaussian()
    fam_gam <- gaussian()
    fam_cgaim <- function(x) x
  } else if(grepl("binom", fam)){
    fam_gam <- fam_gplsim <- binomial()
    fam_gamnl <- fam_binomial(n=1)
  } else if(grepl("poiss", fam)){
    fam_gam <- fam_gplsim <- poisson()
    fam_gamnl <- fam_poisson()
  }
  
  #### gamnl ####
  gamnl_time <- Sys.time()
  fit <- gam_nl(form_gamnl,
                data = dat,
                family = fam_gamnl)
  gamnl_time <- difftime(Sys.time(), gamnl_time, units = "secs")
  
  coefs <- fit$smooth[[1]]$xt$si$B %*% fit$smooth[[1]]$xt$si$alpha ## reparametrise inner coefficients
  
  res <- list(gamnl = list(coefs = c(theta = sign(coefs[1])*coefs/(sqrt(sum(coefs^2))), 
                                     gamma = fit$coefficients[2]),
                           coefs2 = fit$coefficients, 
                           time = gamnl_time,
                           single_test = summary(fit, print = FALSE)$s.table[, "p-value"]))
  
  fit_gam <- gam(y ~ s(ax) + as.factor(Z), 
                 data = dat, family = fam_gam)
  
  res$gam <- list(single_test = summary.gam(fit_gam, print = FALSE)$s.table[, "p-value"])
  yts <- predict(fit, newdata = test, type = "response")
  if(grepl("gauss", fam)){
    res$gamnl["mse"] <- sqrt(sum((test$y-yts[,1])^2)/length(test$y))
  } else if(grepl("binom", fam) || grepl("poiss", fam)){
    res$gamnl["mse"] <- sqrt(sum((test$y-yts)^2)/length(test$y))
  }
  
  ##### gplsim #####
  gplsim_time <- Sys.time()
  fit_gplsim <- gplsim(dat$y,dat$X,dat$Z,user.init=NULL,family = fam_gplsim)       
  gplsim_time <- difftime(Sys.time(), gplsim_time, units = "secs")

  res$gplsim <- list(coefs = c(fit_gplsim$theta, gamma = fit$coefficients[2]),
                     mse = sqrt(sum((fit_gplsim$y-fit_gplsim$fitted.values)^2)/length(fit$y)),
                     time = gplsim_time,
                     single_test = summary(fit_gplsim)$s.table[, "p-value"])
  
  #### cgaim ####
  if(cc == 0)
    return(res)
  if(grepl("binom", fam))
    return(res)
  if(grepl("pois", fam))
    return(res)
  
  dat$yc <- fam_cgaim(dat$y)
  dat$x1 <- dat$X[,1]; dat$x2 <- dat$X[,2]; dat$x3 <- dat$X[,3]
  dat$Z <- as.factor(dat$Z)
  test$yc <- fam_cgaim(test$y)
  test$x1 <- test$X[,1]; test$x2 <- test$X[,2]; test$x3 <- test$X[,3]
  test$Z <- as.factor(test$Z)
  cgaim_time <- Sys.time()
  form_cgaim <- yc ~ g(x1,x2,x3) + Z
  fit <- cgaim(form_cgaim, data = dat)
  cgaim_time <- difftime(Sys.time(), cgaim_time, units = "secs")
  
  yts <- predict(fit, newdata = test)
  
  res$cgaim <- list(coefs = unlist(c(theta = fit$alpha, 
                              gamma = fit$beta[3])),
                           coefs2 = fit$beta, 
                           time = cgaim_time,
                    mse = sqrt(sum((yts-test$y)^2)/length(fit$y)))
  
  return(res)
} ## end fit_models
