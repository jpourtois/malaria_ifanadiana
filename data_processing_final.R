################################################################################################
#### Paper Title: "Climatic, land-use and socio-economic factors can predict malaria dynamics at 
#### fine spatial scales relevant to local health actors: evidence from rural Madagascar"
#### Script title: Data Processing 
#### Script Author: Julie D. Pourtois
#### Updated: July 21st 2022
#### Description: 
################################################################################################

#### 0. Load packages ####

if (!require("dplyr")) install.packages("dplyr");library(dplyr)
if (!require("plyr")) install.packages("plyr");library(plyr)
if (!require("tidyr")) install.packages("tidyr");library(tidyr)
if (!require("stringr")) install.packages("stringr");library(stringr)

#### 1. Set directory ####

setwd("~/Documents/Stanford/Research/Malaria Project/Data/csv_datasets")

#### 2. Load all datasets ####

## 2.1 Fokotany ID and malaria incidence
id_match <- read.csv('fktIDs.csv') # ID to Fokontany match
malaria <- read.csv('dat_izzy_May15.csv') # Main dataset - malaria incidence, population number, distance to health center

## 2.2 Socio-economic variables
bednets <- read.csv('bednets.bs.bf.csv') # Bednet use
wealth <- read.csv('wealth.csv') # Wealth score
road_dist <- read.csv('mdg_fkt_towns_roadDist.csv') # Average distance to road

## 2.3 Climatic variables
alt <- read.csv('altitude.bf.csv') # Altitude
temp_surface <- read.csv('LST.csv') # Surface temperature (GEE)
p <- read.csv('p_tall.csv')

## 2.4 Land use 
land_use <- read.csv('land_use_summary.csv') # Land use in % for 5 classes
deforestation_all <- read.csv('deforestation_all_csv.csv') # Forest cover, gain and loss (2000-2019)
dist_forest <- read.csv('DIST_toForest.csv') # Distance from residential areas to nearest forest fragment
edge <- read.csv('total_class_edges.csv') # Edge length for land-use classes

#### 3. Data processing ####

### 3.0 Malaria incidence 

names(malaria)[names(malaria) == "palu.pc.2nn_0.08"] <- "malaria_total_prop"
names(malaria)[names(malaria) == "palu.und5.pc.3nn_0.19"] <- "malaria_u5_prop"

# Incidence as number of infection per thousand people
malaria$malaria_total_pt <- round(malaria$malaria_total_prop*1000)
malaria$malaria_u5_pt <- round(malaria$malaria_u5_prop*1000)

# Incidence as number of people infected
malaria$malaria_total <- round(malaria$malaria_total_prop*malaria$Population)
malaria$malaria_u5 <- round(malaria$malaria_u5_prop*malaria$Population.u5)

# Time (1-48 months)
malaria$time <- malaria$month + (malaria$year - 2014)*12

# Time labels for figures
malaria <- dplyr::mutate(malaria, 
                  time.labels = paste(month.abb[month], year))

malaria$time.labels <- factor(malaria$time.labels, levels = c("Jan 2014","Feb 2014","Mar 2014","Apr 2014", "May 2014","Jun 2014", "Jul 2014", "Aug 2014", "Sep 2014", "Oct 2014", "Nov 2014", "Dec 2014", "Jan 2015", "Feb 2015",
                                                              "Mar 2015", "Apr 2015", "May 2015", "Jun 2015", "Jul 2015", "Aug 2015", "Sep 2015", "Oct 2015", "Nov 2015", "Dec 2015",
                                                              "Jan 2016", "Feb 2016", "Mar 2016", "Apr 2016", "May 2016", "Jun 2016", "Jul 2016", "Aug 2016", "Sep 2016", "Oct 2016",
                                                              "Nov 2016", "Dec 2016", "Jan 2017", "Feb 2017", "Mar 2017", "Apr 2017", "May 2017", "Jun 2017", "Jul 2017", "Aug 2017",
                                                              "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017"))

### 3.1 Distance to road ##

## Parse commune and Fokontany names and merge with ID file

road_dist <- dplyr::mutate(road_dist,
                    Commune = sapply(str_split(commune_fokotony, "FKT"), "[[", 1),
                    Commune = gsub('.$', '', Commune),
                    Commune = gsub('_', ' ', Commune)) 

road_dist <- dplyr::mutate(road_dist,
                    Fokontany = sapply(str_split(commune_fokotony, "FKT"), "[[", 2),
                    Fokontany = gsub('^.', '', Fokontany),
                    Fokontany = gsub('_', ' ', Fokontany),
                    Fokontany = sapply(str_split(Fokontany, "_"), "[[", 1))
                    
road_dist <- merge(road_dist, id_match, by = c("Fokontany","Commune"))
road_dist <- dplyr::select(road_dist, -c(Fokontany,Commune,commune_fokotony))

### 3.2 Surface temperature

## Parse commune and Fokontany names and merge with ID file
temp_surface$Fokontany <- sub('FKT ', '', temp_surface$LIB_FKT)
temp_surface <- merge(temp_surface, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
temp_surface <- subset(temp_surface, select = -c(LIB_FKT, Fokontany, LIB_COM))

## Create 'suitability index' as squared distance from optimum temperature (see Mordecai et al. 2013)
temp_surface <- dplyr::mutate(temp_surface,
                       LST_C_mean_squared = (LST_C_mean - 25)^2,
                       LST_C_max_squared = (LST_C_max - 25)^2,
                       LST_C_min_squared = (LST_C_min - 25)^2)

### 3.3 Precipitation

## Parse date into year and month
names(p)[names(p) == "OBJEC"] <- "ID"

p <- dplyr::mutate(p,
            year = sapply(strsplit(Date,'_'), "[[", 2),
            year = strtoi(year),
            month = sapply(strsplit(p$Date,'_'), "[[", 3),
            month = as.numeric(month))

p <- subset(p, select = -c(LIB_F,X, Date))

### 3.4 Land use

land_use <- subset(land_use, select = c(ID, X_residenti, X_watermean, X_savannahm, X_ricemean, X_forestmea))
names(land_use) <- c('ID','Residential','Water','Savannah','Rice','Forest')

### 3.5 Deforestation

deforestation_all$Fokontany <- sub('FKT ', '', deforestation_all$LIB_FKT)
deforestation_all <- merge(deforestation_all, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))

## Calculate deforestation over the previous 10 years for each year
deforest_10y <- data_frame('ID' = deforestation_all$ID)

deforest_10y <- dplyr::mutate(deforest_10y,
                       y_2014 = apply(dplyr::select(deforestation_all, X_loss_04me:X_loss_13me), 1, sum, na.rm = TRUE),
                       y_2015 = apply(dplyr::select(deforestation_all, X_loss_05me:X_loss_14me), 1, sum, na.rm = TRUE),
                       y_2016 = apply(dplyr::select(deforestation_all, X_loss_06me:X_loss_15me), 1, sum, na.rm = TRUE),
                       y_2017 = apply(dplyr::select(deforestation_all, X_loss_07me:X_loss_16me), 1, sum, na.rm = TRUE))

deforest_10y <- gather(deforest_10y, year, loss_10y, y_2014:y_2017)
deforest_10y$year <- revalue(deforest_10y$year, c("y_2014"= "2014", "y_2015"="2015", "y_2016" = "2016", 
                                                  "y_2017" = "2017"))
deforest_10y$year <- strtoi(deforest_10y$year)

## Calculate deforestation over the previous 3 years for each year
deforest_3y <- data_frame('ID' = deforestation_all$ID)

deforest_3y <- dplyr::mutate(deforest_3y,
                      y_2014 = apply(dplyr::select(deforestation_all, X_loss_11me:X_loss_13me), 1, sum, na.rm = TRUE),
                      y_2015 = apply(dplyr::select(deforestation_all, X_loss_12me:X_loss_14me), 1, sum, na.rm = TRUE),
                      y_2016 = apply(dplyr::select(deforestation_all, X_loss_13me:X_loss_15me), 1, sum, na.rm = TRUE),
                      y_2017 = apply(dplyr::select(deforestation_all, X_loss_14me:X_loss_16me), 1, sum, na.rm = TRUE))

deforest_3y <- gather(deforest_3y, year, loss_3y, y_2014:y_2017)
deforest_3y$year <- revalue(deforest_3y$year, c("y_2014"= "2014", "y_2015"="2015", "y_2016" = "2016", 
                                                "y_2017" = "2017"))
deforest_3y$year <- strtoi(deforest_3y$year)

### 3.6 Distance from residential area to the forest

dist_forest$Fokontany <- sub('FKT ', '', dist_forest$LIB_FKT)
dist_forest <- merge(dist_forest, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
dist_forest <- subset(dist_forest, select = -c(LIB_FKT, Fokontany, LIB_COM))

### 3.7 Forest edge length

## Parse commune and Fokontany names and merge with ID file
edge <- dplyr::mutate(edge,
               Commune = sapply(str_split(commune_fokotony, "FKT"), "[[", 1),
               Commune = gsub('.$', '', Commune),
               Commune = gsub('_', ' ', Commune),
               Fokontany = sapply(str_split(commune_fokotony, "FKT"), "[[", 2),
               Fokontany = gsub('^.', '', Fokontany),
               Fokontany = gsub('_', ' ', Fokontany),
               Fokontany = sapply(str_split(Fokontany, "_"), "[[", 1))

## Only keep forest edge
edge_forest <- edge[edge$class == 'Forest',]

edge_forest <- merge(id_match, edge_forest, by = c("Fokontany","Commune"), all.x = TRUE)

# Replace NAs with 0 (no data means there was no forest in that Fokontany)
edge_forest$total_edge_m[is.na(edge_forest$total_edge_m)] <- 0

edge_forest <- subset(edge_forest, select = -c(Fokontany,Commune,commune_fokotony, class)) 
names(edge_forest) <- c('ID','edge_forest')

#### 4. Data set merging ####

### 4.1 Merge malaria incidence with Fokontany ID

malaria <- merge(malaria, id_match, by = "ID")

### 4.2 Merge non-temporal datasets (ID)

malaria <- join_all(list(malaria, alt, bednets, land_use, dist_forest, edge_forest, road_dist), by = "ID")

### 4.3 Merge deforestation (ID, year)

malaria <- merge(malaria, deforest_3y, by = c("ID","year"))
malaria <- merge(malaria, deforest_10y, by = c("ID","year"))

### 4.4 Merge wealth (ID, year, month)

malaria <- merge(malaria, wealth[wealth$year != 2018,], by = c("ID","year","month"),all.x = TRUE)

### 4.5 Merge climatic data (ID, year, month)

# First, merge temperature and precipitation
monthly <- merge(temp_surface, p, by = c("ID","year","month"))

# Merge non-lagged temp and precipitation with malaria
malaria <- merge(malaria, monthly, by = c("ID","year","month"))

# Create one-month lag 
monthly_lag <- monthly
monthly_lag$month <- monthly_lag$month + 1
monthly_lag$year[monthly_lag$month == 13] <- monthly_lag$year[monthly_lag$month == 13] + 1
monthly_lag$month[monthly_lag$month == 13] <- 1 
monthly_lag <- rbind(monthly[monthly$year == 2014 & monthly$month == 1,], monthly_lag)
monthly_lag <- monthly_lag[monthly_lag$year != 2018,]
names(monthly_lag) <- c('ID','year','month',"LST_C_mean_lag","LST_C_max_lag","LST_C_min_lag",
                        "LST_C_mean_lag_squared", "LST_C_max_lag_squared","LST_C_min_lag_squared",
                        "Precipitation_lag")

# Create two-month lag
monthly_lag2 <- monthly_lag
monthly_lag2$month <- monthly_lag2$month + 1
monthly_lag2$year[monthly_lag2$month == 13] <- monthly_lag2$year[monthly_lag2$month == 13] + 1
monthly_lag2$month[monthly_lag2$month == 13] <- 1 
monthly_lag2 <- rbind(monthly_lag[monthly_lag$year == 2014 & monthly_lag$month == 1,], monthly_lag2)
monthly_lag2 <- monthly_lag2[monthly_lag2$year != 2018,]
names(monthly_lag2) <- c('ID','year','month',"LST_C_mean_lag2","LST_C_max_lag2","LST_C_min_lag2",
                         "LST_C_mean_lag2_squared", "LST_C_max_lag2_squared", "LST_C_min_lag2_squared",
                         "Precipitation_lag2")

malaria <- merge(malaria, monthly_lag, by = c("ID","year","month"))
malaria <- merge(malaria, monthly_lag2,by = c("ID","year","month"))

## Remove anomalous temperatures
malaria <- subset(malaria, ((!is.na(LST_C_mean_lag)) & (!is.na(LST_C_min_lag)))) 
malaria <- subset(malaria, ((LST_C_mean_lag > 0) & (LST_C_min_lag > -5))) 
malaria$LST_C_max_lag[malaria$LST_C_max_lag > 50] <- 50

malaria <- subset(malaria, ((!is.na(LST_C_mean_lag2)) & (!is.na(LST_C_min_lag2)))) 
malaria <- subset(malaria, ((LST_C_mean_lag2 > 0) & (LST_C_min_lag2 > -5))) 
malaria$LST_C_max_lag2[malaria$LST_C_max_lag2 > 50] <- 50

#### 5. Save data set
write.csv(malaria, 'malaria_v3.csv', row.names=FALSE)

