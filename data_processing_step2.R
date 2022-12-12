### Data processing in preparation of model analysis
# 1. Normalize and scale data
# 2. Average over time for plotting purposes
# 3. Subset to high season for original model and define seasons
# 4. Add additional spatial variables (lat/long, Commune) to assess auto-correlation
# 5. Add malaria lagged by one month for complete data (high and low seasons)

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
  
  # Lat/long and Commune
  commune_info <- read.csv('data/fktIDs.csv')
  coord <- read.csv('data/coordinates.csv')
  
  coord$Fokontany <- sapply(str_split(coord$LIB_FKT, "FKT"), "[[", 2)
  coord$Fokontany <- gsub('^.', '', coord$Fokontany)
  coord$Fokontany <- gsub('_', ' ', coord$Fokontany)
  coord$Fokontany <- sapply(str_split(coord$Fokontany, "_"), "[[", 1)
  
  coord$ID <- commune_info$ID
  
  malaria.norm.high.season <- merge(malaria.norm.high.season, coord[,c('ID','lat','long')], by = 'ID')
  malaria.norm.high.season <- merge(malaria.norm.high.season, commune_info[c('ID', 'Commune')], by = 'ID')
  
  ### 4. Averaged normalized data set with high season only
  high.season.averaged <- aggregate(malaria.norm.high.season,by=list(malaria.norm.high.season$ID), mean, na.rm = TRUE)
  
  ### 5. Total dataset with malaria lag, lat/long and commune
  malaria.norm.total <- malaria.norm
  
  # Add malaria lag
  malaria_pt <- malaria.norm[, c("malaria_total_pt",'time','ID')]
  malaria_pt$time_lag <- malaria_pt$time + 1
  malaria_pt$time <- NULL
  colnames(malaria_pt) <- c('malaria_total_pt_lag','ID','time')
  
  malaria.norm.total <- merge(malaria.norm.total, malaria_pt, by = c('ID','time'))
  malaria.norm.total$malaria_total_pt_lag <- scale(log10(malaria.norm.total$malaria_total_pt_lag + 1))[,1]
  
  # Add lat/long and commune 
  
  malaria.norm.total <- merge(malaria.norm.total, coord[,c('ID','lat','long')], by = 'ID')
  malaria.norm.total <- merge(malaria.norm.total, commune_info[c('ID', 'Commune')], by = 'ID')
  
  
  return(list('averaged.raw' = averaged_data_raw, 'malaria.norm' = malaria.norm, 
              'malaria.norm.high.season' = malaria.norm.high.season, 'high.season.averaged' = high.season.averaged,
              'malaria.norm.total' = malaria.norm.total))
  
}

####  Load and transform data  ####
malaria_raw <- read.csv('output/malaria_v3.csv')

transformed_data <- transformRaw(malaria_raw)

write.csv(transformed_data$averaged.raw, 'output/averaged_raw.csv', row.names = FALSE)
write.csv(transformed_data$malaria.norm, 'output/malaria_norm.csv', row.names = FALSE)
write.csv(transformed_data$malaria.norm.high.season, 'output/malaria_norm_high_season.csv', row.names = FALSE)
write.csv(transformed_data$high.season.averaged, 'output/high_season_averaged.csv', row.names = FALSE)
write.csv(transformed_data$malaria.norm.total, 'output/malaria_norm_total.csv', row.names = FALSE)
