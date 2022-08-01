################################################################################################
#### Paper Title: "Climatic, land-use and socio-economic factors can predict malaria dynamics at 
#### fine spatial scales relevant to local health actors: evidence from rural Madagascar"
#### Script title: Model Selection
#### Script Author: Julie D. Pourtois
#### Updated: July 22nd 2022
#### Description: 
################################################################################################

##### 0. Load packages ####

load_packages <- function() {
  
  if (!require('dplyr')) install.packages('dplyr');library(dplyr)
  if (!require('glmmTMB')) install.packages('glmmTMB');library(glmmTMB)
  if (!require('DHARMa')) install.packages('DHARMa');library(DHARMa)
  if (!require('MuMIn')) install.packages('MuMIn');library(MuMIn)

}

load_packages()

##### 1. Functions #####

transformRaw <- function(malaria_raw) {
  
  ### 1. Average raw data to use for variable maps
  averaged_data_raw <- aggregate(malaria_raw,by=list(malaria_raw$ID), mean, na.rm = TRUE)
  
  ### 2. Normalized data set
  malaria <- malaria_raw
  
  # Normalize data with log function
  malaria <- dplyr::mutate(malaria,
                           Residential = log10(Residential),
                           Rice = log10(Rice),
                           Forest = log10(Forest + 0.001),
                           wscore.n = log10(wscore.n + 0.01),
                           toForest_meanDist = log10(toForest_meanDist + 1),
                           edge_forest = log10(edge_forest + 10),
                           real.dist.csb = log10(real.dist.csb + 10),
                           loss_3y = log10(loss_3y),
                           loss_10y = log10(loss_10y),
                           Precipitation_lag = log10(Precipitation_lag),
                           Precipitation_lag2 = log10(Precipitation_lag2))
  
  # Further scale data
  malaria$real.dist.csb <- scale(malaria[,'real.dist.csb'])
  norm.var <- scale(dplyr::select(malaria, c(alt_bf:Precipitation_lag2)))
  malaria.norm <- cbind(malaria[,c('ID','year','month','time','Population','Population.u5','malaria_total_prop',
                                   'malaria_u5_prop','malaria_total_pt','malaria_u5_pt',
                                   'malaria_total','malaria_u5','real.dist.csb')], norm.var)
  
  malaria.norm <- dplyr::select(malaria.norm, -c('LST_C_mean','LST_C_max','LST_C_min', 
                                            'LST_C_mean_squared','LST_C_max_squared', 'LST_C_min_squared'))
  
  malaria.norm$ID <- as.factor(malaria.norm$ID)
  
  # Only keep rows with no NAs
  malaria.norm <- malaria.norm[complete.cases(malaria.norm),]
  
  ### 3. Normalized data set with high season only
  
  # Top 7 months with high malaria incidence (above 25 cases per thousand people). 
  malaria.norm.high.season <- malaria.norm[(malaria.norm$month > 10 | malaria.norm$month < 6),]
  
  # Assign seasons
  malaria.norm.high.season$season <- NaN
  
  malaria.norm.high.season <- dplyr::mutate(malaria.norm.high.season,
                                    season = replace(season, year == 2014 & month %in% c(1,2,3,4,5), 0),
                                    season = replace(season, year == 2014 & month %in% c(11,12), 1),
                                    season = replace(season, year == 2015 & month %in% c(1,2,3,4,5), 1),
                                    season = replace(season, year == 2015 & month %in% c(11,12), 2),
                                    season = replace(season, year == 2016 & month %in% c(1,2,3,4,5), 2),
                                    season = replace(season, year == 2016 & month %in% c(11,12),3),
                                    season = replace(season, year == 2017 & month %in% c(1,2,3,4,5), 3),
                                    season = replace(season, year == 2017 & month %in% c(11,12), 4))
  
  malaria.norm.high.season$season <- as.factor(malaria.norm.high.season$season)
  
  ### 4. Averaged normalized data set with high season only
  high.season.averaged <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)
  
  return(list('averaged.raw' = averaged_data_raw, 'malaria.norm' = malaria.norm, 
              'malaria.norm.high.season' = malaria.norm.high.season, 'high.season.averaged' = high.season.averaged))
  
}

#### 2. Data loading and final processing ####

setwd("~/Documents/Stanford/Research/Malaria Project/Data/csv_datasets")

malaria_raw <- read.csv('malaria_v3.csv')

transformed_data <- transformRaw(malaria_raw)

averaged.raw <- transformed_data$averaged.raw
malaria.norm <- transformed_data$malaria.norm
malaria.norm.high.season <- transformed_data$malaria.norm.high.season
high.season.averaged <- transformed_data$high.season.averaged

#### 3. Distribution selection ####

## We explore 4 distributions that are appropriate for count data: 
## 1. Poisson distribution
## 2. Zero-inflated Poisson distribution
## 3. Negative-binomial distribution
## 4. Zero-inflated negative-binomial distribution

## 3.1 Poisson distribution
poisson.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared +(1|ID) + (1 |month), 
                       data = malaria.norm.high.season, family = "poisson",na.action="na.exclude")

simulateResiduals(poisson.all, plot = T)
testUniformity(poisson.all) ## KS test and outlier test have p = 0 (Deviation is significant)

## 3.2 ZI Poisson distribution
zipoisson.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared + (1|ID) + (1 |month),
                         data = malaria.norm.high.season, ziformula = ~., family = "poisson",na.action="na.exclude")

simulateResiduals(zipoisson.all, plot = T)
testUniformity(zipoisson.all) ## KS test and outlier test have p = 0 (Deviation is significant)

## 3.3 Negative-binomial distribution
nb.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared +
                         (1|ID) + (1 |month), data = malaria.norm.high.season, family = "nbinom2",na.action="na.exclude")

simulateResiduals(nb.all, plot = T)
testUniformity(nb.all) ## KS test, Dispersion test and outlier test have p = 0 (Deviation is significant)
testZeroInflation(nb.all) # Significant deviation

## 3.4 ZI negative-binomial distribution (ZINB)
zinb.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                    edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                    LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                    Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                    LST_C_min_lag2_squared+ LST_C_max_lag2_squared +
                    (1|ID) + (1 |month), ziformula = ~., data = malaria.norm.high.season, family = "nbinom2",na.action="na.exclude")


simulateResiduals(zinb.all, plot = T)
testUniformity(zinb.all) # No significant deviation
testZeroInflation(zinb.all) # No significant deviation

# ZINB is the clear winner. Let's confirm with AICc.
AICc(poisson.all,zipoisson.all, nb.all, zinb.all)

### --> ZINB will be used for the rest of the analysis.

#### 4. Variable selection ####

### 4.1 First selection 

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

save(glm.model.2, file = 'glm.model.2.v3')
save(aic.matrix.2, file = 'aic.matrix.2.v3')

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

### 4.2  Second step

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

save(glm.model.3, file = 'glm.model.3.v3')
save(aic.matrix.3, file = 'aic.matrix.3.v3')

## Final 'best' model

best.model <- glm.model.3[[which.min(aic.matrix.3)]]

AICc(best.model)
summary(best.model)

save(best.model, file = 'best.model.v3')

## Save model with no random effects for analysis

best.no.random <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb +  
                       wscore.n + loss_3y + Precipitation_lag + LST_C_mean_lag +  
                       LST_C_min_lag + LST_C_mean_lag_squared, 
                     data = malaria.norm.high.season, ziformula = ~., family = "nbinom2",na.action="na.exclude")

save(best.no.random, file = 'best.no.random.v3')

