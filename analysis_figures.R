# Predicting malaria at small spatial scales
# Author: Julie Pourtois
# Co-authors: 


##### Load packages ####

load_packages <- function(){
  
  library('lme4')
  library('plotly')
  library(MASS)
  library(pscl)
  library('lme4')
  library('plotly')
  library('gganimate')
  library('ggplot2')
  library('sf')
  library('tmap')
  library('tmaptools')
  library('leaflet')
  library('transformr')
  library('png')
  library('gifski')
  library('ggsn')
  library('ggpubr')
  library(ape)
  library(lmerTest)
  library('MuMIn')
  library('glmmADMB')
  library(corrr)
  library('car')
  library('glmmTMB')
  library('splines')
  library(piecewiseSEM)
  library(lavaan)
  library(factoextra)
  library(igraph)
  library(tidySEM)
  library(semPlot)
  library(huxtable)
  library(gridExtra)
  library('sf')
  library(DHARMa)
  library(glmulti)
  library(MuMIn)
  library(functClust)
  library(dotwhisker)
  library(ggstatsplot)
  library(parameters)
  library('dplyr')
  library(mapview)
  library(RColorBrewer)
  library(viridis)
  library(tidyr)
  library(raster)
  library(rasterVis)
  library(purrr)
  library(cowplot)
  library(wesanderson)
  library(grid)
  library(gridGraphics)
  library(pscl)
  library('TSdist')
  library('emmeans')
  library(ggeffects)
  
}
load_packages()

##### Functions #####

#Create plot outside of RStudio
myPPlot <- function(thePlot, w=16, h=12, plotPos=c(1,1), invertY=FALSE,
                    newWindow=TRUE, plotDim=c(1,1)) {
  if (newWindow) {
    myWin(w=w, h=h)
    pushViewport(viewport(layout=grid.layout(plotDim[1],plotDim[2])))
  }
  
  plotRow = plotPos[1]
  if (invertY) {
    # invert row order to match MatLab style (y=1 at the bottom, not the top)
    if (plotDim[1]>1) plotRow = plotDim[1] - plotRow + 1
  }
  print(thePlot, vp=viewport(layout.pos.row=plotRow, layout.pos.col=plotPos[2]))
}

#Open new window for plot
myWin <- function(w=16, h=12) {
  quartz(w,h)
  grid.newpage()
}

##################### Data scaling and final processing ######

## Set directory
setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

malaria_raw <- read.csv('malaria.csv')
malaria_raw <- malaria_raw[!is.na(malaria_raw$ID),]

# Average raw data to use for variable maps
averaged_data_raw <- aggregate(malaria_raw,by=list(malaria_raw$ID), mean, na.rm = TRUE)

malaria <- malaria_raw

## Normalize variables

malaria$Residential <- log10(malaria$Residential)
malaria$toRice_meanDist <- log10(malaria$toRice_meanDist)
malaria$Rice <- log10(malaria$Rice)
malaria$Forest <- log10(malaria$Forest + 0.001)
malaria$wscore.n <- log10(malaria$wscore.n + 0.01)
malaria$no_patches_forest <- log10(malaria$no_patches_forest + 0.01)
malaria$forest_patch_area <- log10(malaria$forest_patch_area + 10)
malaria$toForest_meanDist <- log10(malaria$toForest_meanDist + 1)
malaria$toRice_meanDist <- log10(malaria$toRice_meanDist + 0.5)
malaria$edge_forest <- log10(malaria$edge_forest + 10)
malaria$real.dist.csb <- log10(malaria$real.dist.csb + 10)
malaria$loss_3y <- log10(malaria$loss_3y)
malaria$loss_10y <- log10(malaria$loss_10y)
malaria$periToArea <- log10(malaria$periToArea + 0.000000001)
malaria$Precipitation_lag <- log10(malaria$Precipitation_lag)
malaria$Precipitation_lag2 <- log10(malaria$Precipitation_lag2)

## Create scaled dataset ##

malaria$real.dist.csb <- scale(malaria[,'real.dist.csb'])
norm.var <- scale(dplyr::select(malaria, c(alt_bf:wscore.n, periToArea)))
malaria.norm <- cbind(malaria[,c('ID','year','month','time','Population','malaria_total_prop',
                                 'malaria_u5_prop','malaria_total_pt','malaria_u5_pt',
                                 'malaria_total','malaria_u5','real.dist.csb')], norm.var)

malaria.norm$ID <- as.factor(malaria.norm$ID)
malaria.norm <- dplyr::select(malaria.norm, -c('malaria_total_prop_lag2','malaria_u5_prop_lag2','malaria_total_prop_lag','malaria_u5_prop_lag'))

malaria.norm <- malaria.norm[complete.cases(malaria.norm),]


malaria.norm$time.labels <- paste(month.abb[malaria.norm$month],malaria.norm$year)
malaria.norm$time.labels <- factor(malaria.norm$time.labels, levels = c("Apr 2014", "May 2014","Jun 2014", "Jul 2014", "Aug 2014", "Sep 2014", "Oct 2014", "Nov 2014", "Dec 2014", "Jan 2015", "Feb 2015",
                                            "Mar 2015", "Apr 2015", "May 2015", "Jun 2015", "Jul 2015", "Aug 2015", "Sep 2015", "Oct 2015", "Nov 2015", "Dec 2015",
                                            "Jan 2016", "Feb 2016", "Mar 2016", "Apr 2016", "May 2016", "Jun 2016", "Jul 2016", "Aug 2016", "Sep 2016", "Oct 2016",
                                            "Nov 2016", "Dec 2016", "Jan 2017", "Feb 2017", "Mar 2017", "Apr 2017", "May 2017", "Jun 2017", "Jul 2017", "Aug 2017",
                                            "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017"))

cor(dplyr::select(malaria.norm, c(malaria_total_pt, highseason:wscore.n, periToArea)))

## Create averaged dataset ##

# High season only

# Average malaria incidence per month
malaria.per.month <- tapply(malaria.norm$malaria_total_pt, malaria.norm$month, mean, na.rm = TRUE)
#plot(malaria.per.month)

# Top 7 months with high malaria incidence (above 25 cases per thousand people). 
malaria.norm.high.season <- malaria.norm[(malaria.norm$month > 10 | malaria.norm$month < 6),]

# Assign seasons
malaria.norm.high.season$season <- NaN
malaria.norm.high.season$season[malaria.norm.high.season$year == 2014 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 0
malaria.norm.high.season$season[malaria.norm.high.season$year == 2014 & malaria.norm.high.season$month %in% c(11,12)] <- 1
malaria.norm.high.season$season[malaria.norm.high.season$year == 2015 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 1
malaria.norm.high.season$season[malaria.norm.high.season$year == 2015 & malaria.norm.high.season$month %in% c(11,12)] <- 2
malaria.norm.high.season$season[malaria.norm.high.season$year == 2016 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 2
malaria.norm.high.season$season[malaria.norm.high.season$year == 2016 & malaria.norm.high.season$month %in% c(11,12)] <- 3
malaria.norm.high.season$season[malaria.norm.high.season$year == 2017 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 3
malaria.norm.high.season$season[malaria.norm.high.season$year == 2017 & malaria.norm.high.season$month %in% c(11,12)] <- 4

malaria.norm.high.season$season <- as.factor(malaria.norm.high.season$season)

high.season.averaged <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)
high.season.averaged$malaria_total_pt <- scale(log10(high.season.averaged$malaria_total_pt + 0.5))
high.season.averaged$malaria_u5_pt <- scale(log10(high.season.averaged$malaria_u5_pt + 0.5))

##################### Model selection #####################

##### Round 1: Variable selection for lm ####

# All variables to be considered for combinations
variable.vector.1 <- c('toForest_meanDist',
                     'edge_forest','loss_3y','loss_10y', 
                     'LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared', 
                     'LST_C_mean_lag + LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared', 
                     'Precipitation_lag2', 'Precipitation_lag')

# Create matrix of all possible combinations
n <- length(variable.vector.1)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l)
combi.group <- combi.group[-1,]

combi.group <- combi.group[!(combi.group[,7] == 1 & combi.group[,8] == 1),]
combi.group <- combi.group[!(combi.group[,9] == 1 & combi.group[,10] == 1),]

lm.model.1 <- list() # Store list of models
aic.matrix.1 <- rep(0,nrow(combi.group)) # Store vector of AIC

# Evaluate model and store AIC for all combinations
for (row in 1:nrow(combi.group)){
  
  if (row %% 250 == 0) {
    print(paste('Model', row, 'out of', nrow(combi.group)))
    
  }
  
  combi <- combi.group[row,]
  var.select <- variable.vector.1[as.logical(unlist(combi))]
  
  # Variables that are always included in the model
  formula.string <- 'malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + '
  var.string <- paste(var.select, collapse = ' + ')
  formula.string <- paste(formula.string,var.string)
  all.in.formula <- as.formula(formula.string)
  
  lm.model.1[[row]] <- lm(all.in.formula,high.season.averaged)
  aic.matrix.1[row] <- AICc(lm.model.1[[row]])
  
}

hist(aic.matrix.1,20)

# Average top 5% of models
low.quant <- quantile(aic.matrix.1, 0.1)
best.models.1 <- lm.model.1[aic.matrix.1 < low.quant]
avg.best.1 <- model.avg(best.models.1)
summary(avg.best.1)

# Select 'best' model with lowest AICc
best.model.1 <- lm.model.1[[which.min(aic.matrix.1)]]
summary(best.model.1)
vif(best.model.1)
plot(best.model.1)

##### Round 2: Variable selection for glmm, step 1 ####

# All variables to be considered for combinations
variable.vector.2 <- c('toForest_meanDist','edge_forest',
                       'loss_3y','loss_10y',
                       'LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared', 
                       'LST_C_mean_lag + LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared', 
                       'Precipitation_lag2', 'Precipitation_lag')

# Create matrix of all possible combinations
n <- length(variable.vector.2)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l)
combi.group <- combi.group[-1,]

names(combi.group) <- c('forest_dist', 'forest_edge', "loss_3", 'loss_10', 't_lag2', 't_lag1', 'p_lag2', 'p_lag1')

combi.group <- combi.group[!(combi.group[,3] == 1 & combi.group[,4] == 1),] # Loss 3y or 10y
combi.group <- combi.group[combi.group[,5] == 1 | combi.group[,6] == 1,] # Must include temperature
combi.group <- combi.group[!(combi.group[,5] == 1 & combi.group[,6] == 1),] # can't include 1 and 2 month lag
combi.group <- combi.group[combi.group[,7] == 1 | combi.group[,8] == 1,] # Must include precipitation
combi.group <- combi.group[!(combi.group[,7] == 1 & combi.group[,8] == 1),] # can't include 1 and 2 month lag


glm.model.2 <- list() # Store list of models
aic.matrix.2 <- rep(0,nrow(combi.group)) # Store vector of AIC

# Evaluate model and store AIC for all combinations
for (row in 1:nrow(combi.group)){
  
  if (row %% 5 == 0) {
    print(paste('Model', row, 'out of', nrow(combi.group)))
    
  }
  
  combi <- combi.group[row,]
  var.select <- variable.vector.2[as.logical(unlist(combi))]
  
  # Variables that are always included in the model, based on step 1
  formula.string <- 'malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + '
  var.string <- paste(var.select, collapse = ' + ')
  formula.string <- paste(formula.string,var.string)
  formula.string <- paste(formula.string, '+ (1|ID) + (1|month)')
  all.in.formula <- as.formula(formula.string)
  
  glm.model.2[[row]] <- glmmTMB(all.in.formula, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
  aic.matrix.2[row] <- AICc(glm.model.2[[row]])
  
}

save(glm.model.2, file = 'glm.model.2')
save(aic.matrix.2, file = 'aic.matrix.2')

# 'Best' model with lowest AICc

low.quant <- quantile(aic.matrix.2, 0.1)

combi.group.2.best <- combi.group[aic.matrix.2 < low.quant,]
combi.group.2.best$aic <- aic.matrix.2[aic.matrix.2 < low.quant]

best.models.2 <- glm.model.2[aic.matrix.2 < low.quant]
avg.best.2 <- model.avg(best.models.2)
summary(avg.best.2)

best.model.2 <- glm.model.2[[which.min(aic.matrix.2)]]
AICc(best.model.2)
summary(best.model.2)

hist(aic.matrix.2,20)

glm.avg <- model.avg(glm.model.2)
summary(glm.avg)

##### Round 3: Variable selection for glmm, step 2 #####

variable.vector.3 <- c('LST_C_mean_lag','LST_C_min_lag','LST_C_max_lag',
                       'LST_C_mean_lag_squared')

# Create matrix of all possible combinations
n <- length(variable.vector.3)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l)
combi.group <- combi.group[-1,]

combi.group <- combi.group[!(combi.group[,1] == 1 & combi.group[,3] == 1),] # Mean or max
glm.model.3 <- list() # Store list of models
aic.matrix.3 <- rep(0,nrow(combi.group)) # Store vector of AIC

for (row in 1:nrow(combi.group)){
  
  if (row %% 5 == 0) {
    print(paste('Model', row, 'out of', nrow(combi.group)))
    
  }
  
  combi <- combi.group[row,]
  var.select <- variable.vector.3[as.logical(unlist(combi))]
  
  # Variables that are always included in the model, based on step 1
  formula.string <- 'malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + '
  var.string <- paste(var.select, collapse = ' + ')
  formula.string <- paste(formula.string,var.string)
  formula.string <- paste(formula.string, '+ (1|ID) + (1|month)')
  all.in.formula <- as.formula(formula.string)
  
  glm.model.3[[row]] <- glmmTMB(all.in.formula, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
  aic.matrix.3[row] <- AICc(glm.model.3[[row]])
  
}

save(glm.model.3, file = 'glm.model.3')
save(aic.matrix.3, file = 'aic.matrix.3')

## Final 'best' model

best.model <- glm.model.3[[which.min(aic.matrix.3)]]

AICc(best.model)
summary(best.model)

save(best.model, file = 'best.model')

##################### Figures #####################

##### General malaria dynamics #####

sum(malaria$malaria_total, na.rm = TRUE) # Total malaria cases
sum(averaged_data_raw$Population, na.rm = TRUE) # Total population

sum(malaria$malaria_total, na.rm = TRUE)/sum(averaged_data_raw$Population, na.rm = TRUE) # Prop cases per person

# Prop of cases per Fokontany
hist(tapply(malaria$malaria_total_prop, malaria$ID, sum, na.rm = TRUE),20)
range(tapply(malaria$malaria_total_prop, malaria$ID, sum, na.rm = TRUE))

# Prop of cases per Fokontany per month
hist(averaged_data_raw$malaria_total_prop)
mean(averaged_data_raw$malaria_total_prop, na.rm = TRUE)
range(averaged_data_raw$malaria_total_prop, na.rm = TRUE)

# Incidence per month
tapply(malaria$malaria_total_pt, malaria$month, mean, na.rm = TRUE)

mean(averaged_data_raw$malaria_total_pt)
mean(malaria$malaria_total_pt, na.rm = TRUE)

hist(tapply(malaria$malaria_total_prop, malaria$ID, sum, na.rm = TRUE),20)

## Precipitation patterns
tapply(malaria_raw$Precipitation_lag, malaria_raw$month, mean, na.rm = TRUE)

boxplot(Precipitation_lag ~ month, malaria_raw[(malaria_raw$month > 10 | malaria_raw$month < 6) & malaria_raw$year == 2014,])

#Driest month
min(malaria_raw$Precipitation_lag[malaria_raw$month == 11])
max(malaria_raw$Precipitation_lag[malaria_raw$month == 11])
mean(malaria_raw$Precipitation_lag[malaria_raw$month == 11])

# Wettest month
min(malaria_raw$Precipitation_lag[malaria_raw$month == 4])
max(malaria_raw$Precipitation_lag[malaria_raw$month == 4])
mean(malaria_raw$Precipitation_lag[malaria_raw$month == 4])


## Temperature patterns

hist(averaged_data_raw$LST_C_mean_lag)
boxplot(LST_C_mean_lag ~ month, malaria_raw[(malaria_raw$month > 10 | malaria_raw$month < 6)& malaria_raw$year == 2017,])
tapply(malaria_raw$LST_C_mean_lag, malaria_raw$month, mean, na.rm = TRUE)
min(malaria_raw$LST_C_mean_lag[malaria_raw$month == 5])

for (i in 2014:2017){
  df <- malaria_raw[(malaria_raw$month > 10 | malaria_raw$month < 6) & malaria_raw$year == i,]
  print(tapply(df$LST_C_mean_lag, df$month, sd, na.rm = TRUE))
  mean.sd <- mean(tapply(df$LST_C_mean_lag, df$month, sd, na.rm = TRUE))
  
}

hist(tapply(malaria_raw$LST_C_mean_lag[malaria_raw$year == 2016], malaria_raw$ID[malaria_raw$year == 2016], mean, na.rm = TRUE))

##### General analysis #####
load('best.model')

# Poisson dist
poisson.model <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                           wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                           LST_C_min_lag + LST_C_mean_lag_squared + (1|ID) + (1 |month), 
                         data = malaria.norm.high.season, family = "poisson",na.action="na.exclude")

plot(predict(poisson.model, type = 'response'),malaria.norm.high.season$malaria_total_pt)
simulateResiduals(poisson.model, plot = T)
testDispersion(poisson.model)
AICc(poisson.model, best.model)


# Poisson ZI
poissonzi.model <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                           wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                           LST_C_min_lag + LST_C_mean_lag_squared + (1|ID) + (1 |month), 
                         data = malaria.norm.high.season, ziformula = ~.,family = "poisson",na.action="na.exclude")
simulateResiduals(poissonzi.model, plot = T)

AICc(poissonzi.model, best.model)

# NB model
nb.model <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                           wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                           LST_C_min_lag + LST_C_mean_lag_squared + (1|ID) + (1 |month), 
                         data = malaria.norm.high.season, family = "nbinom2",na.action="na.exclude")

plot(predict(nb.model, type = 'response'),malaria.norm.high.season$malaria_total_pt)
simulateResiduals(nb.model, plot = T)
testDispersion(nb.model)

AICc(nb.model, best.model)

# ZINB model

zinb.model <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                      wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                      LST_C_min_lag + LST_C_mean_lag_squared + (1|ID) + (1 |month), 
                    data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")

summary(zinb.model)

plot(predict(best.model, type = 'response'), malaria.norm.high.season$malaria_total_pt)

summary(best.model)
AICc(best.model, zinb.model)
simulateResiduals(best.model, plot = T)
plotResiduals(best.model)
plotQQunif(best.model)

# Model without random effects

no.random <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                       wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                       LST_C_min_lag + LST_C_mean_lag_squared, 
                     data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
summary(no.random)

# Add predictions to dataset

malaria.norm.high.season$pred <- predict(best.model, type = 'response')
malaria.norm.high.season$pred.no.rand <- predict(no.random, type = 'response')

plot(log10(malaria.norm.high.season$malaria_total_pt), log10(malaria.norm.high.season$pred))
plot(log10(malaria.norm.high.season$malaria_total_pt), log10(malaria.norm.high.season$pred.no.rand))

cor(log10(malaria.norm.high.season$malaria_total_pt + 1), log10(malaria.norm.high.season$pred + 1), method = "pearson", use = "complete.obs")
cor(malaria.norm.high.season$malaria_total_pt, malaria.norm.high.season$pred.no.rand, method = "pearson", use = "complete.obs")


# Get RMSE
malaria.norm.high.season$diff <- malaria.norm.high.season$malaria_total_pt - malaria.norm.high.season$pred
rmse.model <- sqrt(mean(malaria.norm.high.season$diff^2))

malaria.norm.high.season$diff.no.rand <- malaria.norm.high.season$malaria_total_pt - malaria.norm.high.season$pred.no.rand
rmse.model.no.rand <- sqrt(mean(malaria.norm.high.season$diff.no.rand^2))


# Calculate rank
high.season.ave <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)
high.season.ave$rank.malaria <- rank(high.season.ave$malaria_total_pt)
high.season.ave$rank.pred <- rank(high.season.ave$pred)
high.season.ave$rank.pred.no.rand <- rank(high.season.ave$pred.no.rand)

top.ID.obs <- high.season.ave$Group.1[rank(high.season.ave$malaria_total_pt) > 175]
top.ID.pred <- high.season.ave$Group.1[rank(high.season.ave$pred) > 175]

top.ID.pred.no.rand <- high.season.ave$Group.1[rank(high.season.ave$pred.no.rand) > 175]

sum(top.ID.obs %in% top.ID.pred)
sum(top.ID.obs %in% top.ID.pred.no.rand)

top.ID.obs <- high.season.ave$Group.1[rank(high.season.ave$malaria_total_pt) > 155]
top.ID.pred <- high.season.ave$Group.1[rank(high.season.ave$pred) > 155]

top.ID.pred.no.rand <- high.season.ave$Group.1[rank(high.season.ave$pred.no.rand) > 155]

sum(top.ID.obs %in% top.ID.pred)
sum(top.ID.obs %in% top.ID.pred.no.rand)


plot(high.season.ave$malaria_total_pt, high.season.ave$pred)
plot(high.season.ave$rank.malaria, high.season.ave$rank.pred.no.rand)

cor(high.season.ave$malaria_total_pt, high.season.ave$pred ,method = "spearman", use = "complete.obs")
cor(high.season.ave$malaria_total_pt, high.season.ave$pred ,method = "pearson", use = "complete.obs")

cor(high.season.ave$malaria_total_pt, high.season.ave$pred.no.rand ,method = "spearman", use = "complete.obs")
cor(high.season.ave$malaria_total_pt, high.season.ave$pred.no.rand ,method = "pearson", use = "complete.obs")

summary(lm(rank(high.season.ave$pred) ~ rank(high.season.ave$malaria_total_pt)))

hist(round(malaria.norm.high.season.pred$pred),40)
hist(malaria.norm.high.season.pred$malaria_total_pt,40)

plot(log10(malaria.norm.high.season.pred$malaria_total_pt + 1),log10(malaria.norm.high.season.pred$pred + 1))


##### Figure 1 - Map of study sites ####

fkt.bound <- st_read('~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT/Limite_FKT_Distr_Ifanadiana.shp')


high.season.pred.averaged <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

high.season.pred.averaged$ID <- as.numeric(as.character(high.season.pred.averaged$Group.1))

map.pred <- inner_join(high.season.pred.averaged, fkt.bound, by = "ID")

write.csv(high.season.pred.averaged[,c('malaria_total_pt','ID')], 'malaria_pred.csv')

centroids <- st_read("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/dd-prediction-pivot-master/data/spatial/study_area/projected-29738/cluster_centroids.shp") %>%
  mutate(cluster_id = as.factor(cluster_id))

roads <- st_read("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/dd-prediction-pivot-master/data/spatial/admin/roads.shp")

#spatial data
commune <- st_read("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/dd-prediction-pivot-master/data/spatial/admin/projected-29738/communes.shp") %>%
  dplyr::select(commune = LIB_COM)

roads <- st_transform(roads, crs(commune)) %>%
  dplyr::select(CLASS, SURFTYPE)
#crop to commune extents
roads <- st_intersection(roads, commune) %>%
  mutate(paved = ifelse(SURFTYPE=="Paved", 1, 0))

road_paved <- roads[roads$paved == 1,]
st_write(road_paved, 'road_paved.shp')

roads_unpaved <- roads[roads$paved == 0,]
st_write(roads_unpaved, 'road_unpaved.shp', append = FALSE)

hc <- st_read("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/dd-prediction-pivot-master/data/spatial/health/Centre-de-santé_Ifanadiana.shp")
hc <- st_transform(hc, crs(commune)) %>%
  dplyr::select(commune = LIB_COM, csb_type = LIB_INFRA, pivot_supp = PIVOT_SUPP)

broad_view <- readPNG('~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/fig1_panelA.png')

land_use_image <- readPNG('~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets/fig1_panelB.png')
plot.new()
rasterImage(land_use_image,0,0,1,1)

fig1A <- ggplot() + background_image(broad_view) + 
  coord_fixed(ratio = 7.16/7.51) +
  theme_minimal()

fig1B <- ggplot() + background_image(land_use_image) + 
  coord_fixed(ratio = 6.88/9.08) +
  theme_minimal()

malaria.per.month <- aggregate(malaria_raw,by=list(malaria_raw$ID,malaria_raw$month), mean, na.rm = TRUE)
climate.per.month <- aggregate(malaria.per.month,by=list(malaria.per.month$month), mean, na.rm = TRUE)
climate.per.month <- dplyr::select(climate.per.month, c(month,LST_C_mean_lag,Precipitation_lag))


malaria.per.month$month <- as.factor(malaria.per.month$month)

fig1C <- ggplot(malaria.per.month, aes(x = month, y = malaria_total_pt)) + 
  geom_boxplot() +
  labs(y = 'Malaria incidence (per thousand)', x = '') +
  scale_x_discrete("Month", breaks = 1:12, labels = month.abb)  +
  theme_minimal()

fig1C

climate.per.month$month <- climate.per.month$month - 1
climate.per.month$month[climate.per.month$month == 0] <- 12

ylim.prim <- c(40, 372)   # in this example, precipitation
ylim.sec <- c(21, 33) # temperature

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1] # there was a bug here

fig1D <- ggplot(climate.per.month, aes(month, Precipitation_lag)) +
  geom_col() +
  geom_line(aes(y = a + LST_C_mean_lag*b), color = "red") +
  scale_y_continuous("Precipitation (mm)", sec.axis = sec_axis(~ (. - a)/b, name = "Temperature (C)")) +
  scale_x_continuous("Month", breaks = 1:12, labels = month.abb)  +
  theme_minimal()


fig1 <- ggarrange(fig1A, fig1B, fig1C, fig1D, labels = c("A","B","C","D"),nrow = 2, ncol = 2, font.label = list(size = 16), vjust = 1)
myPPlot(fig1)



ggsave("fig1.tiff",fig1, units="mm",width=200, height= 150, dpi=300)

##### Figure 2 - Plot of the model coefficients #####
param <- parameters(best.model)

param.both <- param[param$Effects == 'fixed',]
param.both <- param.both[-nrow(param.both),]
param.both <- param.both %>% rename(term = Parameter, estimate = Coefficient, std.error = SE, 
                                    conf.low = CI_low, conf.high = CI_high, 
                                    statistic = z,  p.value = p, model = Component) 

#param.both$estimate[param.both$model == 'zero_inflated'] <- - param.both$estimate[param.both$model == 'zero_inflated']
#param.both$conf.low[param.both$model == 'zero_inflated'] <- - param.both$conf.low[param.both$model == 'zero_inflated']
#param.both$conf.high[param.both$model == 'zero_inflated'] <- - param.both$conf.high[param.both$model == 'zero_inflated']

plot.current <- c("Zero-inflated","Conditional")
names(plot.current) <- c("zero_inflated", "conditional")
fig2<- dwplot(param.both) %>% relabel_predictors(highseason = "Bed net use",
                                          real.dist.csb = "Distance to health center",
                                          wscore.n = "Wealth score",
                                          Residential = "Residential Area",
                                          Rice = "Rice field",
                                          loss_3y = "Forest loss (3 years)",
                                          Precipitation_lag = "Precipitation (1-month lag)",
                                          LST_C_min_lag = "Min LST (1-month lag)",
                                          LST_C_mean_lag = "Mean LST (1-month lag)",
                                          LST_C_mean_lag_squared = "Mean LST index (1-month lag)") + 
  facet_grid(cols = vars(model),scales = "free", labeller = labeller(model = plot.current)) +
  xlab("Estimate (95% CI)") + 
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  theme(legend.position = "none") 

ggsave("fig2.tiff",fig2, units="mm",width=150, height= 80, dpi=300)


ggemmeans(best.model, terms = "LST_C_mean_lag", type = "fe.zi")
fig2.2 <- plot(ggemmeans(best.model, terms = "LST_C_mean_lag", type = "fe.zi")) +
  ggtitle("") +
  xlab("Mean LST (one-month lag, standardized)") +
  ylab("Malaria incidence (per thousand)") +
  theme_minimal()

ggarrange(fig2, fig2.2, widths = c(2, 1.3), labels = c('A','B'))


##### Figure 3 - Plot of the model predictions #####

# Figure 3A
malaria.norm.high.season.pred <- malaria.norm.high.season
malaria.norm.high.season.pred$pred <- predict(best.model, type = 'response')
high.season.pred.averaged <- aggregate(malaria.norm.high.season.pred,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)

fkt.bound <- st_read('~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT/Limite_FKT_Distr_Ifanadiana.shp')

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

high.season.pred.averaged$ID <- as.numeric(as.character(high.season.pred.averaged$Group.1))

map.pred <- inner_join(high.season.pred.averaged, fkt.bound, by = "ID")

plot.pred <- ggplot(map.pred) +
  geom_sf(aes(fill = pred, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,200)) + 
  labs(title = "Predicted incidence", fill = 'Per thousand') +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=10), 
        plot.margin=margin(r=0,l=10, t=20, b=10),
        legend.title=element_text(size=8),legend.text=element_text(size=8))
  
plot.obs <- ggplot(map.pred) +
  geom_sf(aes(fill = malaria_total_pt, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,200)) + 
  labs(title = "Observed incidence", fill = 'Per thousand') +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=10), 
        plot.margin=margin(r=0,l=10, t=20, b=10),
        legend.title=element_text(size=8),legend.text=element_text(size=8))

plot.diff <- ggplot(map.pred) +
  geom_sf(aes(fill = diff, geometry = geometry))+
  blank() +
  labs(title = "Difference", fill = 'Per thousand')  +
  scale_fill_gradient2(limits = c(-35,35)) +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=10), 
        plot.margin=margin(r=0,l=10, t=20, b=10),
        legend.title=element_text(size=8),legend.text=element_text(size=8))

figure3a <- ggarrange(plot.obs, plot.pred, plot.diff, nrow = 1, ncol = 3)
myPPlot(figure3a)

# Figure 3B

fig3b <- ggplot(high.season.ave, aes(x = rank.malaria, y = rank.pred)) +
  geom_point(alpha = 0.7)+
  labs(x = "Malaria incidence (rank)", y = 'Predictions (rank)') +
  theme_minimal()

fig3ab <- ggarrange(figure3a, fig3b, nrow = 1)
fig3ab <- ggarrange(plot.obs, plot.pred, plot.diff, fig3b, labels = c("A","","","B"),nrow = 1, ncol = 4)

# Figure 3C

high.season.pred.time <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$time), mean, na.rm = TRUE)

small.pred <- dplyr::select(malaria.norm.high.season, c('ID','month','year','time.labels','malaria_total_pt','pred'))
long.pred <- gather(small.pred, group, value, c('malaria_total_pt','pred'), factor_key=TRUE)
#long.pred$time <- as.factor(long.pred$time)

CorDistance(high.season.pred.time$malaria_total_prop, high.season.pred.time$pred)

figure3c <- ggplot(long.pred, aes(x=time.labels, y=value, fill=group)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(y = 'Malaria incidence (per thousand)', x = '', fill = '') +
  scale_fill_discrete(labels = c("Observed", "Predicted")) +
  scale_y_continuous(limits = c(0,400)) +
  theme(legend.position="bottom", axis.text.x=element_text(angle = 45, hjust = 1), plot.margin=margin(r=10,l=30, t=20, b=10), axis.title=element_text(size=9))
  

figure1 <- ggarrange(fig3ab, figure3c,labels = c('','C'), nrow = 2, ncol = 1)
myPPlot(figure1)

##### Out of sample predictions #####

high.season.train <- malaria.norm.high.season[malaria.norm.high.season$season %in% c(1,2),]
high.season.test <- malaria.norm.high.season[malaria.norm.high.season$season %in% c(3),]

formula.train <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + 
  wscore.n + loss_3y + LST_C_min_lag + LST_C_mean_lag + LST_C_mean_lag_squared  + Precipitation_lag + (1 | ID) + (1 | month)')

best.model.train <- glmmTMB(formula.train, data = high.season.train, ziformula = ~., family = "nbinom2",na.action="na.exclude")
summary(best.model.train)

high.season.test$pred <- predict(best.model.train, high.season.test, type = 'response')
high.season.test$diff <- high.season.test$malaria_total_pt - high.season.test$pred

rmse.model.test<- sqrt(mean(high.season.test$diff^2))

high.season.test.av <- aggregate(high.season.test,by=list(high.season.test$ID), mean, na.rm = TRUE)

top.ID.obs <- high.season.test.av$Group.1[rank(high.season.test.av$malaria_total_pt) > 175]
top.ID.pred <- high.season.test.av$Group.1[rank(high.season.test.av$pred) > 175]
sum(top.ID.obs %in% top.ID.pred)

top.ID.obs <- high.season.test.av$Group.1[rank(high.season.test.av$malaria_total_pt) > 155]
top.ID.pred <- high.season.test.av$Group.1[rank(high.season.test.av$pred) > 155]
sum(top.ID.obs %in% top.ID.pred)

cor(high.season.test.av$malaria_total_pt,high.season.test.av$pred, method = "spearman")

plot(rank(high.season.test.av$malaria_total_pt),rank(high.season.test.av$pred))

summary(lm(rank(pred) ~ rank(malaria_total_pt), high.season.test.av))

# No random effect

formula.train.fixed <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + 
  wscore.n + loss_3y + LST_C_min_lag + LST_C_mean_lag + LST_C_mean_lag_squared  + Precipitation_lag')

best.model.train.fixed <- glmmTMB(formula.train.fixed, data = high.season.train, ziformula = ~., family = "nbinom2",na.action="na.exclude")
summary(best.model.train.fixed)

high.season.test$pred.fixed <- predict(best.model.train.fixed, high.season.test, type = 'response')
high.season.test$diff.fixed <- high.season.test$malaria_total_pt - high.season.test$pred.fixed

rmse.model.test.fixed <- sqrt(mean(high.season.test$diff.fixed^2))


# Figure

setwd("~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT")
fkt.bound <- st_read('Limite_FKT_Distr_Ifanadiana.shp')
setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

high.season.test.av$ID <- as.numeric(as.character(high.season.test.av$Group.1))

map.pred.test <- inner_join(high.season.test.av, fkt.bound, by = "ID")

plot.pred <- ggplot(map.pred.test) +
  geom_sf(aes(fill = pred, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,262)) + 
  labs(title = "Predicted incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))

plot.obs <- ggplot(map.pred.test) +
  geom_sf(aes(fill = malaria_total_pt, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,262)) + 
  labs(title = "Observed incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))

plot.diff <- ggplot(map.pred.test) +
  geom_sf(aes(fill = diff, geometry = geometry))+
  blank() +
  labs(title = "Difference", fill = 'Per thousand')  +
  scale_fill_gradient2(limits = c(-150,150)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size=12))

figure1c <- ggarrange(plot.obs, plot.pred, plot.diff, nrow = 1, ncol = 3)
myPPlot(figure1c)



##### Figure 4 - SEM #####

psem_model_1 <- psem(
  lm(malaria_total_pt ~ highseason + Rice + Residential + edge_forest + LST_C_mean_lag + 
       Precipitation_lag + real.dist.csb + wscore.n, data = high.season.averaged),
  lm(edge_forest ~ loss_3y, data = high.season.averaged),
  lm(Precipitation_lag ~ alt_bf, data = high.season.averaged),
  lm(LST_C_mean_lag ~ alt_bf, data = high.season.averaged),
  lm(wscore.n ~ roadDist_mean, data = high.season.averaged),
  lm(Residential ~ Savannah, data = high.season.averaged),
  lm(Savannah ~ loss_3y, data = high.season.averaged),
  lm(Rice ~ Precipitation_lag + Savannah, data = high.season.averaged),
  lm(real.dist.csb ~ roadDist_mean, data = high.season.averaged),
  lm(highseason ~ Precipitation_lag + LST_C_mean_lag, data = high.season.averaged),
  lm(loss_3y ~ wscore.n, data = high.season.averaged)
)
p1 <- plot(psem_model_1)
p1

psem_model_1 <- psem(
  lm(malaria_total_pt ~ highseason + Rice + Residential + edge_forest + toForest_meanDist+ LST_C_mean_lag + 
       Precipitation_lag + real.dist.csb + wscore.n + loss_3y, data = high.season.averaged),
  lm(edge_forest ~ loss_3y, data = high.season.averaged),
  lm(Precipitation_lag ~ alt_bf, data = high.season.averaged),
  lm(LST_C_mean_lag ~ alt_bf, data = high.season.averaged),
  lm(wscore.n ~ roadDist_mean, data = high.season.averaged),
  lm(Rice ~ Precipitation_lag, data = high.season.averaged),
  lm(real.dist.csb ~ roadDist_mean, data = high.season.averaged),
  lm(highseason ~ Precipitation_lag + LST_C_mean_lag, data = high.season.averaged),
  lm(loss_3y ~ wscore.n, data = high.season.averaged)
)
psem_model_1



psem_summary <- summary(psem_model_1)

variable_names <- c('Malaria','Bed net','Temperature','Precipitation','Altitude','Forest Edge','Forest loss', ' Wealth score',
                    'Distance to road', 'Distance to HC', 'Residential', 'Distance to forest', ' Rice')

width <- 0.5 + abs(psem_summary$coefficients$Std.Estimate*5)
style <- psem_summary$coefficients$P.Value < 0.05
style[style == TRUE] <- 'normal'
style[style == FALSE] <- 'dashed'

color <- psem_summary$coefficients$Std.Estimate > 0
color[color == TRUE] <- 'LightSeaGreen'
color[color == FALSE] <- 'Salmon2'

label_arrow <- round(psem_summary$coefficients$Std.Estimate,2)

library('DiagrammeR')


predictor_names <- gsub('[.]','',psem_summary$coefficients$Predictor)
response_names <- gsub('[.]','',psem_summary$coefficients$Response)

grViz("

digraph sem_graph {

graph [fontsize = 10, overlap = TRUE]

# several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  malaria_total_pt [label = '@@1-1', width = 5]
  highseason [label = '@@1-2']
  LST_C_mean_lag [label = '@@1-3']
  Precipitation_lag [label = '@@1-4']
  alt_bf [label = '@@1-5']
  edge_forest [label = '@@1-6']
  loss_3y [label = '@@1-7']
  wscoren [label = '@@1-8']
  roadDist_mean [label = '@@1-9']
  realdistcsb [label = '@@1-10']
  Residential [label = '@@1-11']
  toForest_meanDist [label = '@@1-12']
  Rice [label = '@@1-13']
  
  subgraph {
  rank = same; alt_bf; roadDist_mean; toForest_meanDist; Residential;
  }
  
  subgraph {
  rank = same; LST_C_mean_lag; Precipitation_lag; wscoren; realdistcsb;
  }
  
  subgraph {
  rank = same; Rice; highseason;
  }
  
  @@6-1 -> @@7-1 [penwidth = @@2-1, style = @@3-1, color = @@4-1, label = @@5-1, fontcolor = @@4-1]
  @@6-2 -> @@7-2 [penwidth = @@2-2, style = @@3-2, color = @@4-2, label = @@5-2, fontcolor = @@4-2]
  @@6-3 -> @@7-3 [penwidth = @@2-3, style = @@3-3, color = @@4-3, label = @@5-3, fontcolor = @@4-3]
  @@6-4 -> @@7-4 [penwidth = @@2-4, style = @@3-4, color = @@4-4, label = @@5-4, fontcolor = @@4-4]
  @@6-5 -> @@7-5 [penwidth = @@2-5, style = @@3-5, color = @@4-5, label = @@5-5, fontcolor = @@4-5]
  @@6-6 -> @@7-6 [penwidth = @@2-6, style = @@3-6, color = @@4-6, label = @@5-6, fontcolor = @@4-6]
  @@6-7 -> @@7-7 [penwidth = @@2-7, style = @@3-7, color = @@4-7, label = @@5-7, fontcolor = @@4-7]
  @@6-8 -> @@7-8 [penwidth = @@2-8, style = @@3-8, color = @@4-8, label = @@5-8, fontcolor = @@4-8]
  @@6-9 -> @@7-9 [penwidth = @@2-9, style = @@3-9, color = @@4-9, label = @@5-9, fontcolor = @@4-9]
  @@6-10 -> @@7-10 [penwidth = @@2-10, style = @@3-10, color = @@4-10, label = @@5-10, fontcolor = @@4-10]
  @@6-11 -> @@7-11 [penwidth = @@2-11, style = @@3-11, color = @@4-11, label = @@5-11, fontcolor = @@4-11]
  @@6-12 -> @@7-12 [penwidth = @@2-12, style = @@3-12, color = @@4-12, label = @@5-12, fontcolor = @@4-12]
  @@6-13 -> @@7-13 [penwidth = @@2-13, style = @@3-13, color = @@4-13, label = @@5-13, fontcolor = @@4-13]
  @@6-14 -> @@7-14 [penwidth = @@2-14, style = @@3-14, color = @@4-14, label = @@5-14, fontcolor = @@4-14]
  @@6-15 -> @@7-15 [penwidth = @@2-15, style = @@3-15, color = @@4-15, label = @@5-15, fontcolor = @@4-15]
  @@6-16 -> @@7-16 [penwidth = @@2-16, style = @@3-16, color = @@4-16, label = @@5-16, fontcolor = @@4-16]
  @@6-17 -> @@7-17 [penwidth = @@2-17, style = @@3-17, color = @@4-17, label = @@5-17, fontcolor = @@4-17]
  @@6-18 -> @@7-18 [penwidth = @@2-18, style = @@3-18, color = @@4-18, label = @@5-18, fontcolor = @@4-18]
  @@6-19 -> @@7-19 [penwidth = @@2-19, style = @@3-19, color = @@4-19, label = @@5-19, fontcolor = @@4-19]
}

[1]: variable_names
[2]: width
[3]: style
[4]: color
[5]: label_arrow
[6]: predictor_names
[7]: response_names

")

### Compact version



grViz("

digraph sem_graph {

graph [fontsize = 10, overlap = TRUE]

# several 'node' statements
  node [shape = box,
        fontname = Helvetica]
  
  @@1
  
  
}

[1]: psem_summary$psem_summary$coefficients$Predictor
[2]: width
[3]: style
[4]: color
[5]: label_arrow


")

library(data.table)
niv <- variable_names
from <- psem_summary$coefficients$Predictor
to <- psem_summary$coefficients$Response
temp <- data.table(from=factor(from),
                   to=factor(to))

nodes <-   create_node_df(  n=length(niv), label=niv,  width=0.3) 
edges <- create_edge_df(from = temp$from, to = temp$to, 
                        rel = "leading_to")   
graph <- create_graph(  nodes_df = nodes, edges_df = edges)
render_graph(graph)







## Figure S1 - make map of all variables

setwd("~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT")
fkt.bound <- st_read('Limite_FKT_Distr_Ifanadiana.shp')
setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

averaged_data_raw$ID <- as.numeric(as.character(averaged_data_raw$Group.1))

map.raw <- inner_join(averaged_data_raw, fkt.bound, by = "ID")

plot.wealth <- ggplot(map.raw) +
  geom_sf(aes(fill = wscore.n, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Wealth score') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.wealth

plot.dist <- ggplot(map.raw) +
  geom_sf(aes(fill = real.dist.csb, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Dist. to health center') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.dist

plot.bednet <- ggplot(map.raw) +
  geom_sf(aes(fill = highseason, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Bed net use')  +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.bednet

plot.res <- ggplot(map.raw) +
  geom_sf(aes(fill = Residential, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Residential') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.res

plot.rice <- ggplot(map.raw) +
  geom_sf(aes(fill = Rice, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Rice') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.rice

plot.distance.rice <- ggplot(map.raw) +
  geom_sf(aes(fill = toRice_meanDist, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Distance to rice') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.distance.rice

plot.distance.forest <- ggplot(map.raw) +
  geom_sf(aes(fill = toForest_meanDist, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Distance to forest') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.distance.forest

plot.edge.forest <- ggplot(map.raw) +
  geom_sf(aes(fill = edge_forest, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Forest edge') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.edge.forest

plot.peritoarea <- ggplot(map.raw) +
  geom_sf(aes(fill = periToArea, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Perimeter/Area') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.peritoarea

plot.loss.3 <- ggplot(map.raw) +
  geom_sf(aes(fill = loss_3y, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Forest loss (3 years)') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.loss.3

plot.loss.10 <- ggplot(map.raw) +
  geom_sf(aes(fill = loss_10y, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Forest loss (10 years)') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.loss.10

plot.temp.min <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_min_lag2 , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Min temperature') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.min

plot.temp.max <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_max_lag2 , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Max temperature') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.max

plot.temp.mean <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_mean_lag2 , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Mean temperature') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.mean

plot.temp.min.squared <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_min_lag2_squared , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Suitability index (min)') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.min.squared

plot.temp.max.squared <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_max_lag2_squared , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Suitability index (max)') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.max.squared

plot.temp.mean.squared <- ggplot(map.raw) +
  geom_sf(aes(fill = LST_C_mean_lag2_squared , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Suitability index (max)') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.temp.mean.squared

plot.prec <- ggplot(map.raw) +
  geom_sf(aes(fill = Precipitation_lag2 , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4') + 
  labs(fill = 'Precipitation') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))
#plot.prec


myPPlot(ggarrange(plot.wealth, plot.dist, plot.bednet, plot.res, plot.rice, plot.distance.forest,
          plot.distance.rice, plot.edge.forest, plot.peritoarea, plot.loss.10, plot.loss.3, 
          plot.temp.min, plot.temp.max, plot.temp.mean, plot.temp.min.squared, 
          plot.temp.max.squared, plot.temp.min.squared, plot.prec))

## Extra plots - Malaria in different seasons

malaria.0 <- malaria.norm.high.season[malaria.norm.high.season$season == 0,]
malaria.1 <- malaria.norm.high.season[malaria.norm.high.season$season == 1,]
malaria.2 <- malaria.norm.high.season[malaria.norm.high.season$season == 2,]
malaria.3 <- malaria.norm.high.season[malaria.norm.high.season$season == 3,]

averaged.0 <- aggregate(malaria.0,by=list(malaria.0$ID), mean, na.rm = TRUE)
averaged.1 <- aggregate(malaria.1,by=list(malaria.1$ID), mean, na.rm = TRUE)
averaged.2 <- aggregate(malaria.2,by=list(malaria.2$ID), mean, na.rm = TRUE)
averaged.3 <- aggregate(malaria.3,by=list(malaria.3$ID), mean, na.rm = TRUE)

averaged.0$ID <- as.numeric(as.character(averaged.0$Group.1))
averaged.1$ID <- as.numeric(as.character(averaged.1$Group.1))
averaged.2$ID <- as.numeric(as.character(averaged.2$Group.1))
averaged.3$ID <- as.numeric(as.character(averaged.3$Group.1))

top.ID.obs.1 <- averaged.1$Group.1[rank(averaged.1$malaria_total_pt) > 168]
top.ID.pred.1 <- averaged.1$Group.1[rank(averaged.1$pred) > 168]
sum(top.ID.obs.1 %in% top.ID.pred.1)

top.ID.obs.2 <- averaged.2$Group.1[rank(averaged.2$malaria_total_pt) > 176]
top.ID.pred.2 <- averaged.2$Group.1[rank(averaged.2$pred) > 176]
sum(top.ID.obs.2 %in% top.ID.pred.2)

top.ID.obs.3 <- averaged.3$Group.1[rank(averaged.3$malaria_total_pt) > 185]
top.ID.pred.3 <- averaged.3$Group.1[rank(averaged.3$pred) > 185]
sum(top.ID.obs.3 %in% top.ID.pred.3)

map.0 <- inner_join(averaged.0, fkt.bound, by = "ID")
map.1 <- inner_join(averaged.1, fkt.bound, by = "ID")
map.2 <- inner_join(averaged.2, fkt.bound, by = "ID")
map.3 <- inner_join(averaged.3, fkt.bound, by = "ID")

plot.1 <- ggplot(map.1) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = 'Season 1 (2014-2015)', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2 <- ggplot(map.2) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = 'Season 2 (2015-2016)', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.3 <- ggplot(map.3) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = 'Season 3 (2016-2017)', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.1.pred <- ggplot(map.1) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2.pred <- ggplot(map.2) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.3.pred <- ggplot(map.3) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.1.diff <- ggplot(map.1) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2.diff <- ggplot(map.2) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.3.diff <- ggplot(map.3) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

myPPlot(ggarrange(plot.1, plot.2, plot.3, plot.1.pred, plot.2.pred, plot.3.pred, plot.1.diff, plot.2.diff, plot.3.diff, nrow = 3, ncol = 3))






##### Extra analasis ####

p1 <- plot(ggemmeans(best.model, terms = "Rice [n=20]", type = "fixed")) +
  ggtitle("Fixed only")
p2 <- plot(ggemmeans(best.model, terms = "Rice [n=20]", type = "fe.zi")) +
  ggtitle("Fixed  + ZI Effects")
ggarrange(p1,p2)

pred.data <- malaria.norm.high.season[rep(200,50),]
pred.data$real.dist.csb <- seq(min(malaria.norm.high.season$real.dist.csb), max(malaria.norm.high.season$real.dist.csb), length.out = 50)
pred.data$cond_mu <- predict(best.model, newdata = pred.data, type = "response")
pred.data$zi_term <- predict(best.model, type = "zprob", newdata = pred.data)

ggplot(pred.data, aes(x = real.dist.csb)) +
  geom_line(aes(y = cond_mu), color= "red") +
  geom_line(aes(y = zi_term), color = "blue") +
  geom_line(aes(y = cond_mu*(1-zi_term)), color = "green") +
  scale_y_continuous(trans='log10')


simOutZI_soc <- simulateResiduals(best.model, plot = F)
plot(simOutZI_soc)
testZeroInflation(simOutZI_soc)
