load_packages <- function() {
  
  if (!require('tidyr')) install.packages('tidyr');library(tidyr)
  if (!require('dplyr')) install.packages('dplyr');library(dplyr)
  if (!require('glmmTMB')) install.packages('glmmTMB');library(glmmTMB)
  if (!require('huxtable')) install.packages('huxtable');library(huxtable)
  if (!require('sf')) install.packages('sf');library(sf)
  if (!require('png')) install.packages('png');library(png)
  if (!require('ggplot2')) install.packages('ggplot2');library(ggplot2)
  if (!require('ggpubr')) install.packages('ggpubr');library(ggpubr)
  if (!require('grid')) install.packages('grid');library(grid)
  if (!require('ggthemes')) install.packages('ggthemes');library(ggthemes)
  if (!require('parameters')) install.packages('parameters');library(parameters)
  if (!require('dotwhisker')) install.packages('dotwhisker');library(dotwhisker)
  if (!require('emmeans')) install.packages('emmeans');library(emmeans)
  if (!require('ggeffects')) install.packages('ggeffects');library(ggeffects)
  if (!require('piecewiseSEM')) install.packages('piecewiseSEM');library(piecewiseSEM)
  if (!require('DiagrammeR')) install.packages('DiagrammeR');library(DiagrammeR)
  if (!require('DiagrammeRsvg')) install.packages('DiagrammeRsvg');library(DiagrammeRsvg)
  if (!require('magrittr')) install.packages('magrittr');library(magrittr)
  if (!require('rsvg')) install.packages('rsvg');library(rsvg)
  if (!require('MuMIn')) install.packages('MuMIn');library(MuMIn)
  if (!require('DHARMa')) install.packages('DHARMa');library(DHARMa)
  if (!require('stringr')) install.packages('stringr');library(stringr)
  if (!require('ape')) install.packages('ape');library(ape)
  
  
}

load_packages()

averaged.raw <- read.csv('output/averaged_raw.csv')
malaria.norm <- read.csv('output/malaria_norm.csv')
malaria.norm.high.season <- read.csv('output/malaria_norm_high_season.csv')
malaria.norm.total <- read.csv('output/malaria_norm_total.csv')
high.season.averaged <- read.csv('output/high_season_averaged.csv')

#### Assess temporal autocorrelation ####

# Load best model with random effects for Fokontany and month. This model was trained with malaria.norm.high.season
best.model.v3 <- load('output/best.model.v3')
summary(best.model)

# Simulate residuals with DHARMa
sim_res <- simulateResiduals(best.model)
res <- sim_res$fittedResiduals

# Assess temporal auto-correlation for each Fokontany
p_temp <- rep(0,length(unique(malaria.norm.high.season$ID)))
i <- 1

for (ID in unique(malaria.norm.high.season$ID)) {
  data_1 <- malaria.norm.high.season[malaria.norm.high.season$ID == ID,]
  res_1 <- res[malaria.norm.high.season$ID == ID]
  
  p_temp[i] <- testTemporalAutocorrelation(res_1, time = data_1$time)$p.value
  i <- i + 1
}

hist(p_temp,50) # We find strong temporal auto-correlation for most Fokontany.

# Recalculate residuals per time step and assess temporal auto-correlation with DHARMa
sim_res_per_month <- recalculateResiduals(sim_res, group = as.factor(malaria.norm.high.season$time))
testTemporalAutocorrelation(sim_res_per_month, time = unique(malaria.norm.high.season$time)[order(unique(malaria.norm.high.season$time),decreasing = FALSE)])

# It is clear there is some temporal autocorrelation currently not accounted for by the model

#### Assess spatial auto-correlation on original model ####

# By Fokontany
sim_res_per_ID <- recalculateResiduals(sim_res, group = as.factor(malaria.norm.high.season$ID))
testSpatialAutocorrelation(sim_res_per_ID,
                           x = unique(malaria.norm.high.season$long)[order(unique(malaria.norm.high.season$ID),decreasing = FALSE)],
                           y = unique(malaria.norm.high.season$lat)[order(unique(malaria.norm.high.season$ID),decreasing = FALSE)])

# By Commune
sim_res_per_commune <- recalculateResiduals(sim_res, group = as.factor(malaria.norm.high.season$Commune))
testSpatialAutocorrelation(sim_res_per_commune,
                           x = aggregate(long ~ Commune, data = malaria.norm.high.season, FUN = mean)$long,
                           y = aggregate(lat ~ Commune, data = malaria.norm.high.season, FUN = mean)$lat)

# There is significant spatial auto-correlation at the Fokontany level, but not at the commune level

#### Address temporal autocorrelation ####

# Strategy 1: Use Ornstein–Uhlenbeck covariance structure following GLMMTMB vignette
malaria.norm.high.season$time_factor <- numFactor(malaria.norm.high.season$time)

# Note: The following does not converge. We thus attempt to define a OU covariance structure per commune instead. 

# new.model.ou.zi <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
#                              Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
#                              LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + 
#                              (1 | ID) + ou(time_factor-1 |ID), data = malaria.norm.high.season, 
#                            family = "nbinom2", ziformula = ~., na.action = "na.exclude")

new.model.ou.commune <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                                     Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                                     LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + 
                                     (1 | ID) + ou(time_factor-1 |Commune), data = malaria.norm.high.season, 
                                   family = "nbinom2", ziformula = ~., na.action = "na.exclude")

sim_res_per_month <- recalculateResiduals(simulateResiduals(new.model.ou.commune), group = as.factor(malaria.norm.high.season$time))
testTemporalAutocorrelation(sim_res_per_month, time = unique(malaria.norm.high.season$time)[order(unique(malaria.norm.high.season$time),decreasing = FALSE)])

# The covariance structure reduces the temporal autocorrelation but does not address it completely. Next, we directly include 
# malaria lagged by a month. 

# Strategy 2: Directly include lagged malaria. This requires using the whole dataset to avoid losing too many data points. 

new.model.lag <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                       Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                       LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + malaria_total_pt_lag + 
                       (1 | ID) + (1 | month), data = malaria.norm.total, 
                     family = "nbinom2", ziformula = ~., na.action = "na.exclude")

sim_res_per_month <- recalculateResiduals(simulateResiduals(new.model.lag), group = as.factor(malaria.norm.total$time))
testTemporalAutocorrelation(sim_res_per_month, time = unique(malaria.norm.total$time)[order(unique(malaria.norm.total$time),decreasing = FALSE)])

# This does not work either. 

# This means the model with the OU covariance structure is best, even if not perfect. 

#### Address spatial autocorrelation ####

malaria.norm.high.season$pos <- numFactor(scale(malaria.norm.high.season$long), scale(malaria.norm.high.season$lat))
malaria.norm.high.season$group <- 'All'

new.model.ou.commune.matern <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                                  Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                                  LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + 
                                  ou(time_factor-1 |Commune) + mat(pos + 0|group), 
                                data = malaria.norm.high.season, 
                                family = "nbinom2", ziformula = ~., na.action = "na.exclude")

AICc(best.model, new.model.ou.commune.matern) # Decrease in AICc
save(new.model.ou.commune.matern,file = 'output/new.model')


