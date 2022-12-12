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

#### 1. Load data ####

malaria_raw <- read.csv('output/malaria_v3.csv')

averaged.raw <- read.csv('output/averaged_raw.csv')
malaria.norm <- read.csv('output/malaria_norm.csv')
malaria.norm.high.season <- read.csv('output/malaria_norm_high_season.csv')
high.season.averaged <- read.csv('output/high_season_averaged.csv')

#### 2. Distribution selection ####

## We explore 4 distributions that are appropriate for count data: 
## 1. Poisson distribution
## 2. Zero-inflated Poisson distribution
## 3. Negative-binomial distribution
## 4. Zero-inflated negative-binomial distribution

## 2.1 Poisson distribution
poisson.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared +(1|ID) + (1 |month), 
                       data = malaria.norm.high.season, family = "poisson",na.action="na.exclude")

simulateResiduals(poisson.all, plot = T)
testUniformity(poisson.all) ## KS test and outlier test have p = 0 (Deviation is significant)

## 2.2 ZI Poisson distribution
zipoisson.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared + (1|ID) + (1 |month),
                         data = malaria.norm.high.season, ziformula = ~., family = "poisson",na.action="na.exclude")

simulateResiduals(zipoisson.all, plot = T)
testUniformity(zipoisson.all) ## KS test and outlier test have p = 0 (Deviation is significant)

## 2.3 Negative-binomial distribution
nb.all <- glmmTMB(malaria_total_pt ~ highseason + Residential + Rice + real.dist.csb + toForest_meanDist +
                         edge_forest + wscore.n + loss_3y + loss_10y + Precipitation_lag + LST_C_mean_lag +  
                         LST_C_min_lag + LST_C_max_lag + LST_C_mean_lag_squared + LST_C_min_lag_squared+ LST_C_max_lag_squared +
                         Precipitation_lag2 + LST_C_mean_lag2 + LST_C_min_lag2 + LST_C_max_lag2 + LST_C_mean_lag2_squared + 
                         LST_C_min_lag2_squared+ LST_C_max_lag2_squared +
                         (1|ID) + (1 |month), data = malaria.norm.high.season, family = "nbinom2",na.action="na.exclude")

simulateResiduals(nb.all, plot = T)
testUniformity(nb.all) ## KS test, Dispersion test and outlier test have p = 0 (Deviation is significant)
testZeroInflation(nb.all) # Significant deviation

## 2.4 ZI negative-binomial distribution (ZINB)
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

#### 3. Variable selection ####

### 3.1 First selection 

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

### 3.2  Second step

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

