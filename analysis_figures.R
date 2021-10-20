
# Load packages
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
  
}
load_packages()

## Functions

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
myWin <- function(w=16, h=12) {
  quartz(w,h)
  grid.newpage()
}

## Set directory

setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

malaria_raw <- read.csv('malaria.csv')
malaria_raw <- malaria_raw[!is.na(malaria_raw$ID),]

# Average raw data to use for variable maps
averaged_data_raw <- aggregate(malaria_raw,by=list(malaria_raw$ID), mean, na.rm = TRUE)

malaria <- malaria_raw

##################### Data scaling and final processing #####################

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

#cor(dplyr::select(malaria.norm, c(malaria_total_pt, highseason:wscore.n, periToArea)))

## Create averaged dataset ##

# High season only

# Average malaria incidence per month
malaria.per.month <- tapply(malaria.norm$malaria_total_pt, malaria.norm$month, mean, na.rm = TRUE)
plot(malaria.per.month)

# Top 7 months with high malaria incidence (above 25 cases per thousand people). 
malaria.norm.high.season <- malaria.norm[(malaria.norm$month > 10 | malaria.norm$month < 6),]

malaria.norm.high.season$season <- NaN
malaria.norm.high.season$season[malaria.norm.high.season$year == 2014 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 0
malaria.norm.high.season$season[malaria.norm.high.season$year == 2014 & malaria.norm.high.season$month %in% c(11,12)] <- 1
malaria.norm.high.season$season[malaria.norm.high.season$year == 2015 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 1
malaria.norm.high.season$season[malaria.norm.high.season$year == 2015 & malaria.norm.high.season$month %in% c(11,12)] <- 2
malaria.norm.high.season$season[malaria.norm.high.season$year == 2016 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 2
malaria.norm.high.season$season[malaria.norm.high.season$year == 2016 & malaria.norm.high.season$month %in% c(11,12)] <- 3
malaria.norm.high.season$season[malaria.norm.high.season$year == 2017 & malaria.norm.high.season$month %in% c(1,2,3,4,5)] <- 3
malaria.norm.high.season$season[malaria.norm.high.season$year == 2017 & malaria.norm.high.season$month %in% c(11,12)] <- 4

high.season.averaged <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)
high.season.averaged$malaria_total_pt <- scale(log10(high.season.averaged$malaria_total_pt + 0.5))
high.season.averaged$malaria_u5_pt <- scale(log10(high.season.averaged$malaria_u5_pt + 0.5))

##################### Model selection #####################

### Round 1: variable selection

# All variables to be considered for combinations
variable.vector.1 <- c('toRice_meanDist', 'toForest_meanDist',
                     'edge_forest','edge_rice','loss_3y','loss_10y', 'periToArea',
                     'LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2', 
                     'LST_C_mean_lag2_squared + LST_C_min_lag2_squared + LST_C_max_lag2_squared', 
                     'LST_C_mean_lag + LST_C_min_lag + LST_C_max_lag', 
                     'LST_C_mean_lag_squared + LST_C_min_lag_squared + LST_C_max_lag_squared',
                     'Precipitation_lag2', 'Precipitation_lag')

# Create matrix of all possible combinations
n <- length(variable.vector.1)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l)
combi.group <- combi.group[-1,]

lm.model.1 <- list() # Store list of models
aic.matrix.1 <- rep(0,nrow(combi.group)) # Store vector of AIC

# Evaluate model and store AIC for all combinations
for (row in 1:nrow(combi.group)){
  
  if (row %% 2000 == 0) {
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

#hist(aic.matrix.1,20)

# Average top 1% of models
low.quant <- quantile(aic.matrix.1, 0.01)
best.models.1 <- lm.model.1[aic.matrix.1 < low.quant]
avg.best.1 <- model.avg(best.models.1)
summary(avg.best.1)

# Select 'best' model with lowest AICc
best.model.1 <- lm.model.1[[which.min(aic.matrix.1)]]
summary(best.model.1)
vif(best.model.1)
plot(best.model.1)

### Round 2: variable selection

# All variables to be considered for combinations
variable.vector.2 <- c('LST_C_mean_lag', 'LST_C_min_lag','LST_C_max_lag',
                     'LST_C_mean_lag_squared','LST_C_min_lag_squared','LST_C_max_lag_squared')

# Create matrix of all possible combinations
n <- length(variable.vector.2)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l) 
combi.group <- combi.group[-1,]  

lm.model.2 <- list() # Store list of models
aic.matrix.2 <- rep(0,nrow(combi.group)) # Store vector of AIC

# Evaluate model and store AIC for all combinations
for (row in 1:nrow(combi.group)){
  
  combi <- combi.group[row,]
  var.select <- variable.vector.2[as.logical(unlist(combi))]
  
  # Variables that are always included in the model, based on previous step
  formula.string <- 'malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + 
  toForest_meanDist + edge_forest + Precipitation_lag +'
  var.string <- paste(var.select, collapse = ' + ')
  formula.string <- paste(formula.string,var.string)
  all.in.formula <- as.formula(formula.string)
  
  lm.model.2[[row]] <- lm(all.in.formula,high.season.averaged)
  aic.matrix.2[row] <- AICc(lm.model.2[[row]])
  
}

hist(aic.matrix.2,20)

# Select 'best model with lowest AICc
best.model.2 <- lm.model.2[[which.min(aic.matrix.2)]]
summary(best.model.2)
vif(best.model.2)
plot(best.model.2)

### Round 3: Variable selection for glmm

# All variables to be considered for combinations
variable.vector.3 <- c('loss_3y','loss_10y',
                       'LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + LST_C_min_lag2_squared + LST_C_max_lag2_squared', 
                       'LST_C_mean_lag + LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared + LST_C_max_lag_squared', 
                       'Precipitation_lag2', 'Precipitation_lag')

# Create matrix of all possible combinations
n <- length(variable.vector.3)
l <- rep(list(0:1), n)
combi.group <- expand.grid(l)
combi.group <- combi.group[-1,]

combi.group <- combi.group[combi.group[,3] == 1 | combi.group[,4] == 1,] # Must include temperature
combi.group <- combi.group[combi.group[,5] == 1 | combi.group[,6] == 1,] # Must include precipitation


glm.model.3 <- list() # Store list of models
aic.matrix.3 <- rep(0,nrow(combi.group)) # Store vector of AIC

# Evaluate model and store AIC for all combinations
for (row in 1:nrow(combi.group)){
  
  if (row %% 5 == 0) {
    print(paste('Model', row, 'out of', nrow(combi.group)))
    
  }
  
  combi <- combi.group[row,]
  var.select <- variable.vector.3[as.logical(unlist(combi))]
  
  # Variables that are always included in the model, based on step 1
  formula.string <- 'malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + 
  toForest_meanDist + edge_forest + '
  var.string <- paste(var.select, collapse = ' + ')
  formula.string <- paste(formula.string,var.string)
  formula.string <- paste(formula.string, '+ (1|ID) + (1|month) + (1|year)')
  all.in.formula <- as.formula(formula.string)
  
  glm.model.3[[row]] <- glmmTMB(all.in.formula, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
  aic.matrix.3[row] <- AICc(glm.model.3[[row]])
  
}

# 'Best' model with lowest AICc
best.model.3 <- glm.model.3[[which.min(aic.matrix.3)]]
summary(best.model.3)

## Round 4: Variable selection for glmm

# Option 1: 2-month lag for temperature
formula.4.1 <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + 
  toForest_meanDist + edge_forest + loss_3y + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 +
LST_C_mean_lag2_squared + LST_C_min_lag2_squared + LST_C_max_lag2_squared + Precipitation_lag + (1|ID) + (1|month) + (1|year)')

glm.4.1 <- glmmTMB(formula.4.1, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
AICc(glm.4.1) #37141.63

# Option 1: 1-month lag for temperature
formula.4.2 <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + 
  toForest_meanDist + edge_forest + loss_3y + LST_C_mean_lag + LST_C_min_lag + LST_C_max_lag +
LST_C_mean_lag_squared + LST_C_min_lag_squared + LST_C_max_lag_squared + Precipitation_lag + (1|ID) + (1|month) + (1|year)')

glm.4.2 <- glmmTMB(formula.4.2, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
AICc(glm.4.2) #37074.6
summary(glm.4.2)

# Model with lowest AICc and removing mean temperature because of colinearities with max temperature
best.formula.4 <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + wscore.n + 
toForest_meanDist + edge_forest + loss_3y + LST_C_min_lag + LST_C_max_lag +
                            LST_C_min_lag_squared + LST_C_max_lag_squared + Precipitation_lag + (1|ID) + (1|month) + (1|year)')

## Final 'best' model

best.model.4 <- glmmTMB(best.formula.4, data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")
AICc(best.model.4) #37072
summary(best.model.4)

save(best.model.4, file = 'best.model')

##################### Figures #####################

load('best.model')

## Figure 1 - Plot of the model predictions 

best.predictions <- predict(best.model.4, type = 'response')
malaria.norm.high.season.pred <- malaria.norm.high.season
malaria.norm.high.season.pred$pred <- best.predictions
malaria.norm.high.season.pred$diff <- malaria.norm.high.season.pred$malaria_total_pt - malaria.norm.high.season.pred$pred

high.season.pred.averaged <- aggregate(malaria.norm.high.season.pred,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)

setwd("~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT")
fkt.bound <- st_read('Limite_FKT_Distr_Ifanadiana.shp')
setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

high.season.pred.averaged$ID <- as.numeric(as.character(high.season.pred.averaged$Group.1))

map.pred <- inner_join(high.season.pred.averaged, fkt.bound, by = "ID")

plot.pred <- ggplot(map.pred) +
  geom_sf(aes(fill = pred, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,200)) + 
  labs(title = "Predicted incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))
  
plot.obs <- ggplot(map.pred) +
  geom_sf(aes(fill = malaria_total_pt, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,200)) + 
  labs(title = "Observed incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))

plot.diff <- ggplot(map.pred) +
  geom_sf(aes(fill = diff, geometry = geometry))+
  blank() +
  labs(title = "Difference", fill = 'Per thousand')  +
  scale_fill_gradient2(limits = c(-35,35)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size=12))

figure1a <- ggarrange(plot.obs, plot.pred, plot.diff, nrow = 1, ncol = 3)

high.season.pred.time <- aggregate(malaria.norm.high.season.pred,by=list(malaria.norm.high.season$time), mean, na.rm = TRUE)

small.pred <- dplyr::select(malaria.norm.high.season.pred, c('ID','month','year','time','malaria_total_pt','pred'))
long.pred <- gather(small.pred, group, value, c('malaria_total_pt','pred'), factor_key=TRUE)
long.pred$time <- as.factor(long.pred$time)

figure1b <- ggplot(long.pred, aes(x=time, y=value, fill=group)) + 
  geom_boxplot() +
  labs(y = 'Malaria incidence (per thousand)', x = 'Month', fill = '') +
  scale_fill_discrete(labels = c("Observed", "Predicted")) +
  theme(legend.position="bottom")

figure1 <- ggarrange(figure1a, figure1b,labels = c('A','B'),hjust = 0,
                     vjust = 0, nrow = 2, ncol = 1)
myPPlot(figure1)

# Out of sample predictions

high.season.2014.2016 <- malaria.norm.high.season[malaria.norm.high.season$year < 2017,]
high.season.2017 <- malaria.norm.high.season[malaria.norm.high.season$year == 2017,]

formula.2014.2016 <- as.formula('malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + 
  wscore.n + toForest_meanDist + edge_forest + loss_3y + LST_C_min_lag + 
  LST_C_max_lag + LST_C_min_lag_squared + LST_C_max_lag_squared + 
  Precipitation_lag + (1 | ID) + (1 | month)')

best.model.2014.2016 <- glmmTMB(formula.2014.2016, data = high.season.2014.2016, ziformula = ~., family = "nbinom2",na.action="na.exclude")
summary(best.model.2014.2016)

high.season.2017$pred <- predict(best.model.2014.2016, high.season.2017, type = 'response')
high.season.2017$diff <- high.season.2017$malaria_total_pt - high.season.2017$pred

high.season.2017.av <- aggregate(high.season.2017,by=list(high.season.2017$ID), mean, na.rm = TRUE)

setwd("~/Documents/Stanford/Research /Malaria Project/Data/spatial_datasets/Limite FKT")
fkt.bound <- st_read('Limite_FKT_Distr_Ifanadiana.shp')
setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

high.season.2017.av$ID <- as.numeric(as.character(high.season.2017.av$Group.1))

map.pred.2017 <- inner_join(high.season.2017.av, fkt.bound, by = "ID")

plot.pred <- ggplot(map.pred.2017) +
  geom_sf(aes(fill = pred, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,262)) + 
  labs(title = "Predicted incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))

plot.obs <- ggplot(map.pred.2017) +
  geom_sf(aes(fill = malaria_total_pt, geometry = geometry))+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,262)) + 
  labs(title = "Observed incidence", fill = 'Per thousand') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=12))

plot.diff <- ggplot(map.pred.2017) +
  geom_sf(aes(fill = diff, geometry = geometry))+
  blank() +
  labs(title = "Difference", fill = 'Per thousand')  +
  scale_fill_gradient2(limits = c(-150,150)) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size=12))

figure1c <- ggarrange(plot.obs, plot.pred, plot.diff, nrow = 1, ncol = 3)
figure1c

## Figure 2 - Plot of the model coefficients
param <- parameters(best.model.4)

param.both <- param[param$Effects == 'fixed',]
param.both <- param.both[-nrow(param.both),]
param.both <- param.both %>% rename(term = Parameter, estimate = Coefficient, std.error = SE, 
                                            conf.low = CI_low, conf.high = CI_high, 
                                            statistic = z,  p.value = p, model = Component) 

plot.current <- c("Zero-inflated","Conditional")
names(plot.current) <- c("zero_inflated", "conditional")
dwplot(param.both) %>% relabel_predictors(highseason = "Bed net use",
                                              Residential = "Residential Area",
                                              Rice = "Rice field",
                                              real.dist.csb = "Distance to health center",
                                              wscore.n = "Wealth score",
                                              toForest_meanDist = "Distance to forest",
                                              edge_forest = "Forest edge length",
                                              loss_3y = "Forest loss (3 years)",
                                              Precipitation_lag = "Precipitation (1-month lag)",
                                              LST_C_min_lag = "Min LST (1-month lag)",
                                              LST_C_max_lag = "Max LST (1-month lag)",
                                              LST_C_min_lag_squared = "Min LST index (1-month lag)",
                                              LST_C_max_lag_squared = "Max LST index (1-month lag)") + 
  facet_grid(cols = vars(model),scales = "free", labeller = labeller(model = plot.current)) +
  xlab("Estimate (95% CI)") + 
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  theme(legend.position = "none") 


## Figure 3 - SEM

psem_model_1 <- psem(
  lm(malaria_total_pt ~ highseason + Rice + Residential + edge_forest + LST_C_mean_lag + 
       Precipitation_lag + real.dist.csb + wscore.n, data = high.season.averaged),
  lm(edge_forest ~ loss_10y, data = high.season.averaged),
  lm(Precipitation_lag ~ alt_bf, data = high.season.averaged),
  lm(LST_C_mean_lag ~ alt_bf, data = high.season.averaged),
  lm(wscore.n ~ roadDist_mean, data = high.season.averaged),
  lm(Residential ~ Savannah, data = high.season.averaged),
  lm(Savannah ~ loss_10y, data = high.season.averaged),
  lm(Rice ~ Precipitation_lag + Savannah, data = high.season.averaged),
  lm(real.dist.csb ~ roadDist_mean, data = high.season.averaged),
  lm(highseason ~ Precipitation_lag + LST_C_mean_lag, data = high.season.averaged),
  lm(loss_10y ~ wscore.n, data = high.season.averaged)
)
p1 <- plot(psem_model_1)
p1


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

## Extra plots - Malaria in different years

malaria.2014 <- malaria.norm.high.season.pred[malaria.norm.high.season.pred$year == 2014,]
malaria.2015 <- malaria.norm.high.season.pred[malaria.norm.high.season.pred$year == 2015,]
malaria.2016 <- malaria.norm.high.season.pred[malaria.norm.high.season.pred$year == 2016,]
malaria.2017 <- malaria.norm.high.season.pred[malaria.norm.high.season.pred$year == 2017,]

averaged.2014 <- aggregate(malaria.2014,by=list(malaria.2014$ID), mean, na.rm = TRUE)
averaged.2015 <- aggregate(malaria.2015,by=list(malaria.2015$ID), mean, na.rm = TRUE)
averaged.2016 <- aggregate(malaria.2016,by=list(malaria.2016$ID), mean, na.rm = TRUE)
averaged.2017 <- aggregate(malaria.2017,by=list(malaria.2017$ID), mean, na.rm = TRUE)

averaged.2014$ID <- as.numeric(as.character(averaged.2014$Group.1))
averaged.2015$ID <- as.numeric(as.character(averaged.2015$Group.1))
averaged.2016$ID <- as.numeric(as.character(averaged.2016$Group.1))
averaged.2017$ID <- as.numeric(as.character(averaged.2017$Group.1))

map.2014 <- inner_join(averaged.2014, fkt.bound, by = "ID")
map.2015 <- inner_join(averaged.2015, fkt.bound, by = "ID")
map.2016 <- inner_join(averaged.2016, fkt.bound, by = "ID")
map.2017 <- inner_join(averaged.2017, fkt.bound, by = "ID")


plot.2014 <- ggplot(map.2014) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,400)) + 
  labs(title = '2014', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2015 <- ggplot(map.2015) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = '2015', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2016 <- ggplot(map.2016) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = '2016', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2017 <- ggplot(map.2017) +
  geom_sf(aes(fill = malaria_total_pt , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(title = '2017', fill = 'Observed') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7), plot.title = element_text(face = 'bold', hjust = 0.5,size=9)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2014.pred <- ggplot(map.2014) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,400)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2015.pred <- ggplot(map.2015) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2016.pred <- ggplot(map.2016) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2017.pred <- ggplot(map.2017) +
  geom_sf(aes(fill = pred , geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient(low = 'white', high = 'royalblue4',limits = c(0,300)) + 
  labs(fill = 'Predicted') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))




#ggarrange(plot.2014, plot.2015, plot.2016, plot.2017, plot.2014.pred, plot.2015.pred, plot.2016.pred, plot.2017.pred, nrow = 2, ncol = 4)


plot.2014.diff <- ggplot(map.2014) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-250,250)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2015.diff <- ggplot(map.2015) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2016.diff <- ggplot(map.2016) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

plot.2017.diff <- ggplot(map.2017) +
  geom_sf(aes(fill = diff, geometry = geometry), lwd = 0)+
  blank() +
  scale_fill_gradient2(limits = c(-150,150)) +
  labs(fill = 'Difference') +
  theme(legend.title=element_text(size=7),legend.text=element_text(size=7)) + 
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5))

ggarrange(plot.2014, plot.2015, plot.2016, plot.2017, plot.2014.pred, plot.2015.pred, plot.2016.pred, plot.2017.pred, plot.2014.diff, plot.2015.diff, plot.2016.diff, plot.2017.diff, nrow = 3, ncol = 4)

