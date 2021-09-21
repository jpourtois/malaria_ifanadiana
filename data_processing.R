# Load packages
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

### Read main data sets ###

# Set directory

setwd("~/Documents/Stanford/Research /Malaria Project/Data/csv_datasets")

## Read datasets and prepare for merging ##

id_match <- read.csv('fktIDs.csv') # ID to Fokontany match

malaria <- read.csv('dat_izzy_May15.csv') # Main dataset - malaria incidence, population number, distance to health center
malaria <- merge(malaria, id_match, by = "ID")

#coordinates <- read.csv('coordinates.csv') # Lat-long and other Fokontany coordinates
#coordinates$Fokontany <- sub('FKT ', '', coordinates$LIB_FKT)
#coordinates <- merge(coordinates, id_match,by= "Fokontany")
#coordinates <- subset(coordinates, select = -c(LIB_FKT, Fokontany, Commune, x_unknown, y_unknown))
#coordinates$ID[duplicated(coordinates$ID)]

alt <- read.csv('altitude.bf.csv') # Altitude
bednets <- read.csv('bednets.bs.bf.csv') # Bednet use

wealth <- read.csv('wealth.csv')

# Land use dataset

land_use <- read.csv('land_use_summary.csv')
land_use <- subset(land_use, select = c(ID, X_residenti, X_watermean, X_savannahm, X_ricemean, X_forestmea))
names(land_use) <- c('ID','Residential','Water','Savannah','Rice','Forest')

# Protected forest

protected_f <- read.csv('MDG_FKT_remainingForest_protected.csv')
protected_f$Commune <- sapply(str_split(protected_f$commune_fokotony, "FKT"), "[[", 1)
protected_f$Commune <- gsub('.$', '', protected_f$Commune)
protected_f$Commune <- gsub('_', ' ', protected_f$Commune)
protected_f$Fokontany <- sapply(str_split(protected_f$commune_fokotony, "FKT"), "[[", 2)
protected_f$Fokontany <- gsub('^.', '', protected_f$Fokontany)
protected_f$Fokontany <- gsub('_', ' ', protected_f$Fokontany)
protected_f$Fokontany <- sapply(str_split(protected_f$Fokontany, "_"), "[[", 1)
protected_f <- merge(protected_f, id_match, by = c("Fokontany","Commune"))
protected_f <- subset(protected_f, select = -c(Fokontany,Commune,commune_fokotony)) 
names(protected_f)[names(protected_f) == "sum"] <- "protec_sum"
# Deforestation

deforestation_all <- read.csv('deforestation_all_csv.csv') # Forest cover, gain and loss (2000-2019)
deforestation_all$Fokontany <- sub('FKT ', '', deforestation_all$LIB_FKT)
deforestation_all <- merge(deforestation_all, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))

# Deforestation over the previous 10 years
deforest_10y <- data_frame('ID' = deforestation_all$ID)
deforest_10y$y_2014 <- apply(select(deforestation_all, X_loss_04me:X_loss_13me), 1, sum, na.rm = TRUE)
deforest_10y$y_2015 <- apply(select(deforestation_all, X_loss_05me:X_loss_14me), 1, sum, na.rm = TRUE)
deforest_10y$y_2016 <- apply(select(deforestation_all, X_loss_06me:X_loss_15me), 1, sum, na.rm = TRUE)
deforest_10y$y_2017 <- apply(select(deforestation_all, X_loss_07me:X_loss_16me), 1, sum, na.rm = TRUE)
deforest_10y <- gather(deforest_10y, year, loss_10y, y_2014:y_2017)
deforest_10y$year <- revalue(deforest_10y$year, c("y_2014"= "2014", "y_2015"="2015", "y_2016" = "2016", 
                                                    "y_2017" = "2017"))
deforest_10y$year <- strtoi(deforest_10y$year)

# Deforestation over the previous three years
deforest_3y <- data_frame('ID' = deforestation_all$ID)
deforest_3y$y_2014 <- apply(select(deforestation_all, X_loss_11me:X_loss_13me), 1, sum, na.rm = TRUE)
deforest_3y$y_2015 <- apply(select(deforestation_all, X_loss_12me:X_loss_14me), 1, sum, na.rm = TRUE)
deforest_3y$y_2016 <- apply(select(deforestation_all, X_loss_13me:X_loss_15me), 1, sum, na.rm = TRUE)
deforest_3y$y_2017 <- apply(select(deforestation_all, X_loss_14me:X_loss_16me), 1, sum, na.rm = TRUE)
deforest_3y <- gather(deforest_3y, year, loss_3y, y_2014:y_2017)
deforest_3y$year <- revalue(deforest_3y$year, c("y_2014"= "2014", "y_2015"="2015", "y_2016" = "2016", 
                                                            "y_2017" = "2017"))
deforest_3y$year <- strtoi(deforest_3y$year)


# Deforestation over the previous year

deforestation <- subset(deforestation_all, select = c(ID, X_loss_13me, X_loss_14me, X_loss_15me, X_loss_16me, X_loss_17me))
deforestation <- gather(deforestation, year, loss, X_loss_13me:X_loss_17me)
deforestation$year <- revalue(deforestation$year, c("X_loss_13me"="2013", "X_loss_14me"="2014", "X_loss_15me" = "2015", 
                                                    "X_loss_16me" = "2016","X_loss_17me" = "2017"))
deforestation$year <- strtoi(deforestation$year)
loss_of <- deforestation[deforestation$year != 2013,]
loss_previous <- deforestation[deforestation$year != 2017,]
loss_previous$year <- loss_previous$year + 1
names(loss_of)[names(loss_of) == "loss"] <- "loss_of"
names(loss_previous)[names(loss_previous) == "loss"] <- "loss_previous"

# Surface temperature
temp_surface <- read.csv('LST.csv') # Surface temperature (GEE)
temp_surface$Fokontany <- sub('FKT ', '', temp_surface$LIB_FKT)
temp_surface <- merge(temp_surface, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
temp_surface <- subset(temp_surface, select = -c(LIB_FKT, Fokontany, LIB_COM))

temp_surface$LST_C_mean_squared <- (temp_surface$LST_C_mean - 25)^2
temp_surface$LST_C_max_squared <- (temp_surface$LST_C_max - 25)^2
temp_surface$LST_C_min_squared <- (temp_surface$LST_C_min - 25)^2

# NDWI
ndwi <- read.csv('NDWI.csv') # Soil moisture
ndwi$Fokontany <- sub('FKT ', '', ndwi$LIB_FKT)
ndwi <- merge(ndwi, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
ndwi <- subset(ndwi, select = -c(LIB_FKT, Fokontany, LIB_COM))

# Distance between residential areas and rice field
dist_rice <- read.csv('DIST_toRice.csv') # Distance from residential areas to rice fields
dist_rice$Fokontany <- sub('FKT ', '', dist_rice$LIB_FKT)
dist_rice <- merge(dist_rice, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
dist_rice <- subset(dist_rice, select = -c(LIB_FKT, Fokontany, LIB_COM))

# Distance between residential areas and forest
dist_forest <- read.csv('DIST_toForest.csv') # Distance from residential areas to rice fields
dist_forest$Fokontany <- sub('FKT ', '', dist_forest$LIB_FKT)
dist_forest <- merge(dist_forest, id_match,by.x = c("Fokontany", "LIB_COM"),by.y = c("Fokontany", "Commune"))
dist_forest <- subset(dist_forest, select = -c(LIB_FKT, Fokontany, LIB_COM))

# Forest edge

edge <- read.csv('total_class_edges.csv')
edge$Commune <- sapply(str_split(edge$commune_fokotony, "FKT"), "[[", 1)
edge$Commune <- gsub('.$', '', edge$Commune)
edge$Commune <- gsub('_', ' ', edge$Commune)
edge$Fokontany <- sapply(str_split(edge$commune_fokotony, "FKT"), "[[", 2)
edge$Fokontany <- gsub('^.', '', edge$Fokontany)
edge$Fokontany <- gsub('_', ' ', edge$Fokontany)
edge$Fokontany <- sapply(str_split(edge$Fokontany, "_"), "[[", 1)
edge <- merge(edge, id_match, by = c("Fokontany","Commune"))
edge <- subset(edge, select = -c(Fokontany,Commune,commune_fokotony)) 
edge_wide <- spread(edge, key = 'class', value = 'total_edge_m')
names(edge_wide)[names(edge_wide) == "Forest"] <- "edge_forest"
names(edge_wide)[names(edge_wide) == "Residential area"] <- "edge_res"
names(edge_wide)[names(edge_wide) == "Rice field"] <- "edge_rice"
names(edge_wide)[names(edge_wide) == "Savanna"] <- "edge_savanna"
edge_wide$`Water bodies` <- NULL

# Forest fragmentation
forest_frag <- read.csv('forestFrag_fkt_all.csv')
forest_frag$Commune <- sapply(str_split(forest_frag$commune_fokotony, "FKT"), "[[", 1)
forest_frag$Commune <- gsub('.$', '', forest_frag$Commune)
forest_frag$Commune <- gsub('_', ' ', forest_frag$Commune)
forest_frag$Fokontany <- sapply(str_split(forest_frag$commune_fokotony, "FKT"), "[[", 2)
forest_frag$Fokontany <- sapply(str_split(forest_frag$Fokontany, ".tif"), "[[", 1)
forest_frag$Fokontany <- gsub('^.', '', forest_frag$Fokontany)
forest_frag$Fokontany <- gsub('_', ' ', forest_frag$Fokontany)
forest_frag$Fokontany <- sapply(str_split(forest_frag$Fokontany, "_"), "[[", 1)
forest_frag <- merge(forest_frag, id_match, by = c("Fokontany","Commune"))
forest_frag <- subset(forest_frag, select = -c(Fokontany,Commune,commune_fokotony)) 

# Distance to road

road_dist <- read.csv('mdg_fkt_towns_roadDist.csv')
road_dist$Commune <- sapply(str_split(road_dist$commune_fokotony, "FKT"), "[[", 1)
road_dist$Commune <- gsub('.$', '', road_dist$Commune)
road_dist$Commune <- gsub('_', ' ', road_dist$Commune)
road_dist$Fokontany <- sapply(str_split(road_dist$commune_fokotony, "FKT"), "[[", 2)
road_dist$Fokontany <- gsub('^.', '', road_dist$Fokontany)
road_dist$Fokontany <- gsub('_', ' ', road_dist$Fokontany)
road_dist$Fokontany <- sapply(str_split(road_dist$Fokontany, "_"), "[[", 1)
road_dist <- merge(road_dist, id_match, by = c("Fokontany","Commune"))
road_dist <- select(road_dist, -c(Fokontany,Commune,commune_fokotony))

p <- read.csv('p_tall.csv')
names(p)[names(p) == "OBJEC"] <- "ID"
p$year <- sapply(strsplit(p$Date,'_'), "[[", 2)
p$month <- sapply(strsplit(p$Date,'_'), "[[", 3)
p$year <- strtoi(p$year)
p$month <- as.numeric(p$month) 
p <- subset(p, select = -c(LIB_F,X, Date))

## Merge non-temporal datasets ##

malaria <- merge(malaria, alt, by = "ID")
malaria <- merge(malaria, bednets, by = "ID")
malaria <- merge(malaria, land_use, by = "ID")
malaria <- merge(malaria, protected_f, by = "ID")
malaria <- merge(malaria, dist_rice, by = "ID")
malaria <- merge(malaria, dist_forest, by = "ID")
malaria <- merge(malaria, forest_frag, by = "ID")
malaria <- merge(malaria, edge_wide, by = "ID")
malaria <- merge(malaria, road_dist, by = "ID")

## Merge deforestation ##

malaria <- merge(malaria, loss_of, by = c("ID","year"))
malaria <- merge(malaria, loss_previous, by = c("ID","year"))
malaria <- merge(malaria, deforest_3y, by = c("ID","year"))
malaria <- merge(malaria, deforest_10y, by = c("ID","year"))

## Merge monthly climate/weather data ##

# First, merge all monthly data we want

monthly <- merge(ndwi,temp_surface, by = c("ID","year","month"))
monthly <- merge(monthly, p, by = c("ID","year","month"))
monthly <- merge(monthly, malaria[,c("palu.pc.2nn_0.08","palu.und5.pc.3nn_0.19","ID","year","month")],by = c("ID","year","month"))


# Create lag for all monthly

monthly_lag <- monthly
monthly_lag$month <- monthly_lag$month + 1
monthly_lag$year[monthly_lag$month == 13] <- monthly_lag$year[monthly_lag$month == 13] + 1
monthly_lag$month[monthly_lag$month == 13] <- rep(1,length(monthly_lag$month[monthly_lag$month == 13]))
monthly_lag <- rbind(monthly[monthly$year == 2014 & monthly$month == 1,], monthly_lag)
monthly_lag <- monthly_lag[monthly_lag$year != 2018,]
names(monthly_lag) <- c('ID','year','month',"NDWI_mean_lag","NDWI_max_lag","NDWI_min_lag","LST_C_mean_lag",
                        "LST_C_max_lag","LST_C_min_lag","LST_C_mean_lag_squared", "LST_C_max_lag_squared",
                        "LST_C_min_lag_squared","Precipitation_lag","malaria_total_prop_lag","malaria_u5_prop_lag")

monthly_lag2 <- monthly_lag
monthly_lag2$month <- monthly_lag2$month + 1
monthly_lag2$year[monthly_lag2$month == 13] <- monthly_lag2$year[monthly_lag2$month == 13] + 1
monthly_lag2$month[monthly_lag2$month == 13] <- rep(1,length(monthly_lag2$month[monthly_lag2$month == 13]))
monthly_lag2 <- rbind(monthly_lag[monthly_lag$year == 2014 & monthly_lag$month == 1,], monthly_lag2)
monthly_lag2 <- monthly_lag2[monthly_lag2$year != 2018,]
names(monthly_lag2) <- c('ID','year','month',"NDWI_mean_lag2","NDWI_max_lag2","NDWI_min_lag2","LST_C_mean_lag2",
                         "LST_C_max_lag2","LST_C_min_lag2","LST_C_mean_lag2_squared", "LST_C_max_lag2_squared",
                         "LST_C_min_lag2_squared","Precipitation_lag2","malaria_total_prop_lag2","malaria_u5_prop_lag2")

malaria <- merge(malaria, monthly_lag, by = c("ID","year","month"))
malaria <- merge(malaria, monthly_lag2,by = c("ID","year","month"))

# Merge wealth
malaria <- merge(malaria, wealth[wealth$year != 2018,], by = c("ID","year","month"),all.x = TRUE)

## Additional data processing ##

# Incidence as proportion of population
names(malaria)[names(malaria) == "palu.pc.2nn_0.08"] <- "malaria_total_prop"
names(malaria)[names(malaria) == "palu.und5.pc.3nn_0.19"] <- "malaria_u5_prop"
names(malaria)[names(malaria) == "mean_forestPatch_area_sqm"] <- "forest_patch_area"

# Incidence as number of infection per thousand people
malaria$malaria_total_pt <- round(malaria$malaria_total_prop*1000)
malaria$malaria_u5_pt <- round(malaria$malaria_u5_prop*1000)

malaria$malaria_total_pt_lag <- round(malaria$malaria_total_prop_lag*1000)
malaria$malaria_u5_pt_lag <- round(malaria$malaria_u5_prop_lag*1000)

# Incidence as number of people infected
malaria$malaria_total <- round(malaria$malaria_total_prop*malaria$Population)
malaria$malaria_u5 <- round(malaria$malaria_u5_prop*malaria$Population.u5)

# Time
malaria$time <- malaria$month + (malaria$year - 2014)*12

## Remove messed up temperatures

malaria <- malaria[malaria$LST_C_mean_lag2 > 0,]
malaria <- malaria[malaria$LST_C_min_lag2 > -5,]
#malaria <- malaria[malaria$LST_C_max_lag2 < 50,]
malaria$LST_C_max_lag2[malaria$LST_C_max_lag2 > 50] <- 50

## Replace NAs with 0 for some forest variables

malaria$no_patches_forest[is.na(malaria$no_patches_forest)] <- 0
malaria$forest_patch_area[is.na(malaria$forest_patch_area)] <- 0
malaria$total_forest_area[is.na(malaria$total_forest_area)] <- 0
malaria$edge_forest[is.na(malaria$edge_forest)] <- 0

## Calculate forest perimeter to area ratio

malaria$periToArea <- malaria$edge_forest/malaria$total_forest_area
malaria$periToArea[malaria$periToArea == Inf] <- 0
malaria$periToArea[is.na(malaria$periToArea)] <- 0

write.csv(malaria, 'malaria.csv')
