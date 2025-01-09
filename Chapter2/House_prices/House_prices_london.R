rm(list = ls())

library(gamFactory)
library(mgcViz)
library(gamlss)
library(scoringRules)
library(gridExtra)
library(viridis)
load("./data/house_prices_data_london.rdata")

dat <- dat[dat$city == "LONDON",]
test <- dat[dat$date>"2022-11-30",]
dat <- dat[dat$date<="2022-11-30",]

## define formula 
formulas <- list(std = list(log_price ~ s(grid_e, grid_n, k = 12^2) +
                                  s(ah3ahah) + s(imd) + s(log_dist_tube) +
                                   property_type + old_new + duration + category_type
                                 ,~ property_type + category_type + s(imd)
                                 ,~1,~1)) # skewness and kurtosis are constant

formulas <- c(formulas, 
              lapply(c(3,5,7,8,10,12)^2, 
                     function(ii) list(formula(paste("log_price ~ s_nest(ks, m = c(4, 2), trans = trans_mgks())",
                                         "s_nest(ahah_measures, m = c(4, 2), trans = trans_linear()) ",
                                         "s(imd) "," s(log_dist_tube)",paste(" s(grid_e,grid_n, k =", ii,") "),
                                         "property_type "," old_new "," duration "," category_type", sep = "+"))
                                       ,~ property_type + category_type + s(imd)
                                       ,~1,~1))
              )

names(formulas) <- c("std", "s(e,n)","s(ah)", paste0("k", c(3,5,7,8,10,12)))

# !!! mclapply create a 3 copies of dat
if (file.exists("fitted_house_london.rda")) {
  load("fitted_house_london_measures.rda")
} else {
  start <- format(Sys.time())
  fits <- parallel::mclapply(formulas,
                             function(ff) {
                               start <- format(Sys.time())
                               fit <- gam_nl(ff,family = fam_shash(),
                                      data = dat, optimizer="efs")
                               time <- difftime(format(Sys.time()),start,  units = "min")
                               return(list(fit = fit, time = time))
                               },mc.cores = min(4, length(formulas)))
  fitting_time <- lapply(fits, function(ff) ff$time)
  fits <- lapply(fits, function(ff) ff$fit)
  fits <- lapply(fits, getViz)
  names(fits) <- names(formulas)
  save(fits, fitting_time, file = "fitted_house_london_measures.rda")
}

summary(fits[[1]])
summary(fits[[2]])

perf <- data.frame(
  mod = names(formulas),
  logloss = -sapply(fits, function(xx) logLik(xx)),
  mse = sapply(fits, function(xx) sum((predict(xx, newdata = test)[,1]-test$log_price)^2)),
  R2 = sapply(fits, FUN = function(ff) 1-var(residuals(ff, type = "response"))/var(ff$y)),
  AIC = sapply(fits, AIC),
  df = sapply(lapply(fits, "[[", "edf1"), sum),
  time = sapply(fitting_time, as.numeric, units = "mins"))

library(tidyverse)
th <- theme(axis.text = element_text(size=15),
            legend.title = element_text(size=17),
            axis.title=element_blank(),
            legend.text = element_text(size=15),
            title = element_text(size=17))
newDF <- perf %>% 
  pivot_longer(col= !contains(c("mod", "logloss", "df")), 
               values_to = "column_values", 
               names_to = "column_names")
newDF$column_names <- factor(newDF$column_names, levels = unique(newDF$column_names))
newDF %>% ggplot(aes(x = factor(mod, levels = perf$mod), y = column_values)) +
  geom_point() +
  facet_wrap(vars(column_names), scales = "free_y", nrow = 1)+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        axis.title = element_blank())+th


## Outer and inner terms' plots
xy <- lapply(1:length(fits[[4]]$smooth), 
       function(ii){
          sm <- lapply(fits[-(1:3)], function(ff) {
            ss <- sm(ff, ii)
            outer <- data.frame(x = plot(ss)$data$fit$x,
                               y = plot(ss)$data$fit$y,
                               kk = sqrt(ff$smooth[[5]]$df+1))
          if(class(ss)[1]%in%c("mgks", "si")){
            dd <- plot(ss, inner=TRUE)$data$fit
            inner <- data.frame(x = dd$x,
                               y = dd$y*sign(dd$y[1]),
                               se = dd$se,
                               kk = sqrt(ff$smooth[[5]]$df+1))
            outer$x <- outer$x*sign(dd$y[1])
          } else {inner <- NULL}
          return(list(outer = outer, 
                      inner = inner))  
          }) # end sm_lapply
          list(outer= do.call("rbind", lapply(sm, "[[", "outer")),
               inner = do.call("rbind", lapply(sm, "[[", "inner")))         
       })# end xy lapply
names(xy) <- sapply(fits[[4]]$smooth, "[[", "label")
names(xy)[c(1,4,5)] <- c("s(mgks)", "s(d)", "s(lon,lat)")

# plot of outer effects
p <- do.call("grid.arrange", 
        c(lapply(names(xy)[!(names(xy) == "s(lon,lat)")],
                 function(ss){
  ggplot(xy[[ss]]$outer, aes(x=x, y=y, color=factor(kk))) + geom_line()+
                     ggtitle(ss) + scale_color_viridis(discrete = TRUE) +
                     scale_fill_viridis(discrete = TRUE)+ labs(color='k')+
                     theme_light()
}), ncol=3))

# plot of inner coefficients
pd <- position_dodge(1)
p1 <- lapply(c("s(mgks)", "s(ahah_measures)"),
             function(ss){
               ggplot(xy[[ss]]$inner, aes(x=factor(kk),y=y, color=factor(kk))) + 
                 geom_errorbar(aes(ymax = y + se, ymin = y - se), position = pd)+
                 geom_point(position=pd)+
                 ggtitle(ss) + scale_color_viridis(discrete = TRUE) +
                 facet_wrap(facets = vars(x),scales = "free_y", ncol=5)+
                 scale_fill_viridis(discrete = TRUE)+ labs(color='Dim')+
                 theme_light()+xlab("k")
             })

library(wacolors)
terms <- predict(fits[[6]], type="terms")[, -(1:6)]
colnames(terms) <- names(xy)
p3 <- do.call("grid.arrange", 
        c(lapply(colnames(terms),
                 function(ss){
                   ggplot(data.frame(x=dat$grid_e, y=dat$grid_n, z=terms[, ss]), 
                          aes(x=x, y=y, color=z)) + geom_point(size = 0.5)+
                     ggtitle(ss)+ scale_color_wa_c("olympic", midpoint = 0)+theme_light()
                 }), ncol=3))
fit <- fits[[6]]

#### single-effect plots ####
get_plot <- function(ns){
  th <- theme(axis.text = element_text(size=22),
              axis.title=element_text(size=25),
              legend.title = element_text(size=22),
              legend.text = element_text(size=18),
              title = element_text(size=22))
  plot(fit, select = ns) + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciPoly() + l_fitLine() + th
}

get_plot(1) + xlab(expression(tilde(s)[mgks](y[-i]))) + ylab(expression(s(tilde(s)[mgks](y[-i]))))
get_plot(2) + xlab(expression(ahah^T~alpha[ahah])) + ylab(expression(s(ahah^T~alpha[ahah])))
plot(fit, select = 2, inner=TRUE)  +
  theme(axis.text = element_text(size=22),
        axis.title=element_text(size=25),
        legend.title = element_text(size=22),
        legend.text = element_text(size=18),
        title = element_text(size=22))+
  geom_abline(linetype = 2, slope = 0, intercept = 0, alpha=0.5)+
  ylab(expression(alpha[ahah])) + xlab("") +
  scale_x_discrete(labels=c("Retail", "Health", "Physical env.", "Air quality")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5))
get_plot(3) + xlab("imd") + ylab(expression(s(imd)))
get_plot(4) + xlab("log_distance") + ylab(expression(s(d)))

plot(fit, select = 6, all.Terms=TRUE) + scale_x_discrete(limits = levels(dat$property_type))
plot(fit, select = 7, all.Terms=TRUE)
plot(fit, select = 8, all.Terms=TRUE)
plot(fit, select = 9, all.Terms=TRUE)
plot(fit, select = 10, allTerms=TRUE)

### check 2d
check2D(fits[[2]], type = "y", x1 = dat$grid_e, x2 = dat$grid_n, maxpo = Inf) +
  l_gridCheck2D(gridFun = median, stand = FALSE) +  l_points()
check2D(fits[[2]], type = "response", x1 = dat$grid_e, x2 = dat$grid_n, maxpo = Inf) +
  l_gridCheck2D(gridFun = median, stand = FALSE, bw = 2000) + l_points()

