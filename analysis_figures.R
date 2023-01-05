################################################################################################
#### Paper Title: "Climatic, land-use and socio-economic factors can predict malaria dynamics at 
#### fine spatial scales relevant to local health actors: evidence from rural Madagascar"
#### Script title: Analysis and figures 
#### Script Author: Julie D. Pourtois
#### Updated: August 1st 2022
#### Description: 
################################################################################################

#### O. LOAD PACKAGES ####

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
  if (!require('lmtest')) install.packages('lmtest');library(lmtest)
  
}

load_packages()

#### 1. FUNCTIONS ####

### 1.1 Analysis

analysis_variables <- function(df.raw, df.averaged.raw) {
  
  ### 1. Malaria incidence patterns
  
  total_malaria <- sum(df.raw$malaria_total, na.rm = TRUE) # Total malaria cases
  total_population <- sum(df.averaged.raw$Population, na.rm = TRUE) # Total population
  
  average_per_capita <- total_malaria/total_population # Prop cases per person
  
  # Prop. of cases per Fokontany
  range_average_fok <- range(tapply(df.raw$malaria_total_prop, df.raw$ID, sum, na.rm = TRUE))
  
  # Prop of cases per Fokontany per month
  mean_monthly_fok <- mean(df.averaged.raw$malaria_total_pt, na.rm = TRUE)
  range_monthly_fok <- range(df.averaged.raw$malaria_total_pt, na.rm = TRUE)
  
  # Number of missing data
  
  missing_row <- 0
  
  for (id in unique(df.raw$ID)) {
    
    missing_row <- missing_row + (48 - nrow(df.raw[df.raw$ID == id,]))
    
  }
  
  missing_row_ave <- missing_row/length(unique(df.raw$ID))
  
  malaria_analysis_table <- hux(Metric = c('Total malaria', 'Total population','Average infections PC','Min. infections P.C.','Max infections P.C.',
                                           'Mean monthly P.T.','Min. monthly P.T.', 'Max. monthly P.T','Mean missing data'), 
                                Number = c(total_malaria, total_population, average_per_capita, range_average_fok[1], range_average_fok[2],
                                           mean_monthly_fok, range_monthly_fok[1], range_monthly_fok[2], missing_row_ave))  %>% 
    set_bold(1, everywhere)
  
  malaria_analysis_table
  
  # Incidence per month
  monthly_malaria <- tapply(df.raw$malaria_total_pt, df.raw$month, mean, na.rm = TRUE)
  
  ### 2. Precipitation and temperature patterns
  
  ## Precipitation patterns
  p_monthly_mean <- tapply(df.raw$Precipitation, df.raw$month, mean, na.rm = TRUE)

  ## Temperature patterns
  temp_monthly_mean <- tapply(df.raw$LST_C_mean[df.raw$LST_C_mean > 0], df.raw$month[df.raw$LST_C_mean > 0], mean, na.rm = TRUE)
  
  malaria_climate_per_month <- hux('Month' = c('January','February','March','April','May','June','July','August','September','October','November','December'), 
                                   'Malaria incidence (p.t.)' = monthly_malaria,
                                   'Precipitation (mm)' = p_monthly_mean,
                                   'Mean temperature (C)' = temp_monthly_mean) %>% 
    set_bold(1, everywhere)
  
  return(list('malaria' = malaria_analysis_table, 'climate' = malaria_climate_per_month))
  
  
}

in_sample_pred <- function(model, df.norm.high.season) {
  
  ### 1. Within-sample predictions
  
  df.norm.high.season$pred <- predict(model, type = 'response')
 
  # Correlation between predicted and observed
  cor.pred <- cor(log10(df.norm.high.season$malaria_total_pt + 1), log10(df.norm.high.season$pred + 1), method = "pearson", use = "complete.obs")
  
  # Get RMSE
  df.norm.high.season$diff <- df.norm.high.season$malaria_total_pt - df.norm.high.season$pred
  rmse.model <- sqrt(mean(df.norm.high.season$diff^2))
  
  # Average predictions over time
  
  high.season.ave <- aggregate(df.norm.high.season,by=list(df.norm.high.season$ID), mean, na.rm = TRUE)
  
  cor.spatial.spearman <- cor(high.season.ave$malaria_total_pt, high.season.ave$pred ,method = "spearman", use = "complete.obs")
  cor.spatial.pearson <- cor(high.season.ave$malaria_total_pt, high.season.ave$pred ,method = "pearson", use = "complete.obs")
  
  # Rank according to malaria incidence and predictions
  high.season.ave$rank.malaria <- rank(high.season.ave$malaria_total_pt)
  high.season.ave$rank.pred <- rank(high.season.ave$pred)
 
  # Rank prediction accuracy for top 10% Fokontany for malaria incidence
  top.ID.obs <- high.season.ave$Group.1[rank(high.season.ave$malaria_total_pt) > 175]
  top.ID.pred <- high.season.ave$Group.1[rank(high.season.ave$pred) > 175]
  
  top_20 <- sum(top.ID.obs %in% top.ID.pred)
  
  # Rank prediction accuracy for top 20% Fokontany for malaria incidence
  top.ID.obs.40 <- high.season.ave$Group.1[rank(high.season.ave$malaria_total_pt) > 155]
  top.ID.pred.40 <- high.season.ave$Group.1[rank(high.season.ave$pred) > 155]
  
  top_40 <- sum(top.ID.obs.40 %in% top.ID.pred.40)
  
  ### Make table 
 
  model_analysis <- hux(Metrics = c("RMSE","Pearson's R", 
                                    "Pearson's R, spatial", "Spearman's rho, spatial",
                                    "Top 20 (10%) rank", "Top 40 (20%) rank"),
                        Results = c(rmse.model, cor.pred, 
                                    cor.spatial.pearson, cor.spatial.spearman,
                                    top_20, top_40))
  model_analysis
  
}

out_sample_pred <- function(model, df.norm.high.season){
  
  ### Out-of sample predictions
  
  # Define training and testing data
  high.season.train <- df.norm.high.season[df.norm.high.season$season %in% c(1,2),]
  high.season.test <- df.norm.high.season[df.norm.high.season$season %in% c(3),]
  
  model.formula <- model$call$formula
  model.train <- glmmTMB(model.formula, data = high.season.train, ziformula = ~., family = "nbinom2",na.action="na.exclude")
  
  # Predict on test data
  high.season.test$pred <- predict(model.train, high.season.test, type = 'response')
  high.season.test$diff <- high.season.test$malaria_total_pt - high.season.test$pred
  
  # Calculate RMSE
  rmse.model.test<- sqrt(mean(high.season.test$diff^2))
  
  # Compare RMSE to mean and range of malaria incidence
  high.season.malaria.mean <- mean(df.norm.high.season$malaria_total_pt)
  high.season.malaria.range <- range(df.norm.high.season$malaria_total_pt)
  
  # Calculate spatial rank
  high.season.test.av <- aggregate(high.season.test,by=list(high.season.test$ID), mean, na.rm = TRUE)
  high.season.test.av$rank.malaria <- rank(high.season.test.av$malaria_total_pt)
  high.season.test.av$rank.pred <- rank(high.season.test.av$pred)
  
  # Rank prediction accuracy for top 10% Fokontany for malaria incidence
  top.ID.obs.out.20 <- high.season.test.av$Group.1[rank(high.season.test.av$malaria_total_pt) > 175]
  top.ID.pred.out.20 <- high.season.test.av$Group.1[rank(high.season.test.av$pred) > 175]
  top.20.out <- sum(top.ID.obs.out.20 %in% top.ID.pred.out.20)
  
  # Rank prediction accuracy for top 20% Fokontany for malaria incidence
  top.ID.obs.out.40 <- high.season.test.av$Group.1[rank(high.season.test.av$malaria_total_pt) > 155]
  top.ID.pred.out.40 <- high.season.test.av$Group.1[rank(high.season.test.av$pred) > 155]
  top.40.out <- sum(top.ID.obs.out.40 %in% top.ID.pred.out.40)
  
  spearman.out <- cor(high.season.test.av$malaria_total_pt,high.season.test.av$pred, method = "spearman")
  pearson.out <- cor(high.season.test.av$malaria_total_pt,high.season.test.av$pred, method = "pearson")
  
  
  ### Make table 
  
  model_analysis <- hux(Metrics = c("RMSE, out-of-sample",
                                    "High-season malaria average", "High-season malaria min.", "High-season malaria max.",
                                    "Top 20 (10%) rank, out-of-sample", "Top 40 (20%) rank, out-of-sample",
                                    "Pearson's R, out-of-sample", "Spearman's rho, out-of-sample"),
                        Results = c(rmse.model.test,
                                    high.season.malaria.mean, high.season.malaria.range[1], high.season.malaria.range[2],
                                    top.20.out, top.40.out,
                                    pearson.out, spearman.out))
  model_analysis
  
}

calculateMoran <- function(model, df.norm.high.season){
  
  # Get residuals with model predictions
  df.norm.high.season$pred <- predict(model, type = 'response')
  df.norm.high.season$res <- df.norm.high.season$malaria_total_pt - df.norm.high.season$pred
  
  # Aggregate residuals over time to get sum of residuals for each Fokontany
  res.by.time <- aggregate(res ~ ID, data = df.norm.high.season, FUN = sum, na.rm = TRUE)
  
  # Get lat and long for each Fokontany
  lats <- aggregate(lat~ID, data = df.norm.high.season, FUN = mean, na.rm = TRUE)
  long <- aggregate(long~ID, data = df.norm.high.season, FUN = mean, na.rm = TRUE)
  
  # Calculate distance matrix for Moran's I test
  fok.dist <- as.matrix(dist(cbind(long$long, lats$lat)))
  
  fok.dist.inv <- 1/fok.dist
  diag(fok.dist.inv) <- 0
  
  moran.test <- Moran.I(res.by.time$res, fok.dist.inv)
  return(moran.test)
  
}

calculateACF <- function(model, df.norm.high.season){
  
  # Get residuals from model predictions
  df.norm.high.season$pred <- predict(model, type = 'response')
  df.norm.high.season$res <- df.norm.high.season$malaria_total_pt - df.norm.high.season$pred
  
  # Aggregate residuals over space to get sum of residuals for each time point
  res.by.time <- aggregate(res ~ time, data = df.norm.high.season, FUN = sum, na.rm = TRUE)
  
  # Calculate ACF
  acf.test <- acf(res.by.time$res)
  
  # Test for autocorrelation with DW test
  dw.test <- lmtest::dwtest(res.by.time$res ~ 1, order.by = res.by.time$time)
  
  return(list(acf = acf.test, test = dw.test))
  
}

### 1.2 Tables 

table_model_coefficients <- function(model){
  
  param <- parameters(model)
  
  param.both <- param[param$Effects == 'fixed',]
  param.both$df_error <- NULL
  param.both$Group <- NULL
  param.both$Effects <- NULL
  param.both <- param.both %>% dplyr::rename(term = Parameter, estimate = Coefficient, std.error = SE, 
                                             conf.low = CI_low, conf.high = CI_high, 
                                             statistic = z,  p.value = p, model = Component) 
  
  param.cond <- param.both[param.both$model == 'conditional',]
  param.cond$variable <- c('Intercept', 'Bed net use', 'Residential area (log10)', 'Rice field area (log10)',
                       'Distance to h.c. (log10)', 'Wealth score (log10)', 'Forest loss (log10)', 
                       'Precipitation, 1-month lag (log10)', 'Mean LST, 1-month lag', 'Min LST, 1-month lag',
                       'Mean LST index, 1-month lag')
  
  param.zi <- param.both[param.both$model == 'zero_inflated',]
  param.zi$variable <- c('Intercept', 'Bed net use', 'Residential area (log10)', 'Rice field area (log10)',
                         'Distance to h.c. (log10)', 'Wealth score (log10)', 'Forest loss (log10)', 
                         'Precipitation, 1-month lag (log10)', 'Mean LST, 1-month lag', 'Min LST, 1-month lag',
                         'Mean LST index, 1-month lag')
  
  coef.table <- hux(
    Variable = param.cond$variable,
    cond.est  = exp(param.cond$estimate),
    cond.ci =  paste(round(exp(param.cond$conf.low), digits = 2),'-', round(exp(param.cond$conf.high), digits =2)),
    zi.est = exp(param.zi$estimate),
    zi.ci = paste(round(exp(param.zi$conf.low), digits = 2),'-', round(exp(param.zi$conf.high), digits =2))
  )
  
  coef.table <- coef.table %>% 
    set_contents(1, 2:5, c("Estimate", "CI (95%)", "Estimate", "CI (95%)")) %>% 
    insert_row("", "Conditional","", "Zero-inflated", "", after = 0) %>% 
    merge_cells(1, 2:3) %>% 
    merge_cells(1, 4:5) %>% 
    set_bold(1, everywhere) %>% 
    set_bold(2, everywhere) %>%
    set_align(1, everywhere, "center") %>%
    set_align(2, everywhere, "center") %>%
    set_bottom_border(row = 2, col = everywhere) 
  
  coef.table
  
  quick_docx(coef.table, file = 'output/coef_table.docx')
  #quick_pdf(coef.table, file = 'coef_table.pdf')
  
}

table_model_avg <- function(){
  
  # Load all models
  load('output/glm.model.2.v3') 
  load('output/aic.matrix.2.v3')
  
  # Average of top 10% models after step 2
  low.quant <- quantile(aic.matrix.2, 0.1)
  
  best.models.2 <- glm.model.2[aic.matrix.2 < low.quant]
  avg.best.2 <- model.avg(best.models.2)
  summary(avg.best.2)
  
  top10_model <- parameters(avg.best.2)
  
  average_model <- hux(Variable = top10_model$Parameter,
      Estimate = top10_model$Coefficient,
      SE = top10_model$SE,
      'p-value' = top10_model$p)
  
  average_model <- average_model %>% 
    set_bold(1, everywhere)
  
  quick_docx(average_model, file = 'output/average_model.docx')
  
}

### 1.3 Figures

fig_map_malaria <- function(df.averaged.raw){
  
  ## Data processing for qGIS

  # Export observation data
  write.csv(df.averaged.raw[,c('malaria_total_pt','ID')], 'output/malaria_obs_v3.csv')
  
  # Load Fokontany borders
  fkt.bound <- st_read('data/Limite_FKT/Limite_FKT_Distr_Ifanadiana.shp')
  fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))
  
  # Load road data
  roads <- st_read('data/roads/roads.shp')
  
  # Load commune border data
  commune <- st_read('data/commune_projected_29738/communes.shp') %>%
    dplyr::select(commune = LIB_COM)
  
  # Load road data
  roads <- st_transform(roads, st_crs(commune)) %>%
    dplyr::select(CLASS, SURFTYPE)
  
  #crop to commune extents
  roads <- st_intersection(roads, commune) %>%
    mutate(paved = ifelse(SURFTYPE=="Paved", 1, 0))
  
  # Select and export paved roads
  road_paved <- roads[roads$paved == 1,]
  #st_write(road_paved, 'output/road_paved_v3.shp')
  
  # Select and export unpaved roads
  roads_unpaved <- roads[roads$paved == 0,]
  #st_write(roads_unpaved, 'output/road_unpaved_v3.shp', append = FALSE)
  
  ## Make Figure 1A and 1B
  
  # Load Figure 1A (made in qGIS)
  broad_view <- readPNG('output/fig1A_v3.png')
  
  # Load Figure 1B (made in qGIS)
  land_use_image <- readPNG('output/fig1_panelBalt.png')
  
  fig1A <- ggplot() + background_image(broad_view) + 
    coord_fixed(ratio = 17.27/21.34) + # Change image ratio if needed
    theme_minimal()
  
  fig1B <- ggplot() + background_image(land_use_image) + 
    coord_fixed(ratio = 16.2/20.44) + # Change image ratio if needed
    theme_minimal()
  
  fig1AB <- ggarrange(fig1A, fig1B, labels = c("A","B"), font.label = list(size = 10))
  
  # Save tiff file
  ggsave("figures/fig1.pdf",fig1AB, units="mm",width=225, height= 112, dpi=300)
  ggsave("figures/Fig1.tiff",fig1AB, units="in",width=5.2, height= 2.4, dpi=300)
  
}

fig_coefficients <- function(model, df.norm.high.season, df.raw){
  
  ## Coefficients for conditional model
  
  param <- parameters(model)
  
  param.both <- param[param$Effects == 'fixed',]
  param.both <- param.both[-nrow(param.both),]
  param.both$df_error <- NULL
  param.both$Group <- NULL
  param.both$Effects <- NULL
  param.both <- param.both %>% dplyr::rename(term = Parameter, estimate = Coefficient, std.error = SE, 
                                             conf.low = CI_low, conf.high = CI_high, 
                                             statistic = z,  p.value = p, model = Component) 
  
  param.cond <- param.both[param.both$model == 'conditional',]
  plot.current <- c("Conditional")
  names(plot.current) <- c("conditional")
  
  param.cond$estimate <- exp(param.cond$estimate)
  param.cond$CI <- exp(param.cond$estimate)
  param.cond$conf.low <- exp(param.cond$conf.low) 
  param.cond$conf.high <- exp(param.cond$conf.high)
  
  fig2a <- dwplot(param.cond) %>% relabel_predictors(highseason = "Bed net use",
                                                     real.dist.csb = "Distance to health center (log10)",
                                                     wscore.n = "Wealth score (log10)",
                                                     Residential = "Residential area (log10)",
                                                     Rice = "Rice field area (log10)",
                                                     loss_3y = "Three-year forest loss (log10)",
                                                     Precipitation_lag = "Precipitation, 1-month lag (log10)",
                                                     LST_C_min_lag = "Min LST, 1-month lag",
                                                     LST_C_mean_lag = "Mean LST, 1-month lag",
                                                     LST_C_mean_lag_squared = "Mean LST index, 1-month lag") + 
    facet_grid(cols = vars(model),scales = "free", labeller = labeller(model = plot.current)) +
    xlab("Exp. of Estimate (95% CI)") + 
    geom_vline(xintercept = 1, colour = "grey60", linetype = 2) + 
    theme_bw() +
    theme(legend.position = "none") 
  
  ### Figure 2B
  
  ## Calculate marginal effect for LST by hand (to consider both LST and LST_squared), credit to Michelle Evans
  
  ID.eg <- 156
  time.eg <- 37
  
  # Create data frame with mean of each variable
  pred.temp.data <- df.norm.high.season %>%
    dplyr::select(-malaria_total_pt, -ID, -month, -season, - Commune, - time_factor, -pos, -group, -LST_C_mean_lag, -LST_C_mean_lag_squared) %>%
    summarise_all(mean) %>%
    mutate(ID = ID.eg, time = time.eg)
  
  # Re-add data that can't be averaged
  extra.data <- df.norm.high.season[df.norm.high.season$ID == ID.eg & df.norm.high.season$time == time.eg,
                                  c('season', 'Commune','time_factor','pos','group')]
  
  pred.temp.data <- cbind(pred.temp.data, extra.data)
  
  # Create data frame with just LST and LST index
  temp.join <- df.norm.high.season %>%
    dplyr::select(LST_C_mean_lag, LST_C_mean_lag_squared) %>%
    mutate(LST_C_mean_lag = round(LST_C_mean_lag,1),
           LST_C_mean_lag_squared = round(LST_C_mean_lag_squared,1)) %>%
    distinct() %>%
    arrange(LST_C_mean_lag) %>% 
    group_by(LST_C_mean_lag) %>%
    arrange(LST_C_mean_lag_squared) %>%
    slice(1)
  
  # Replicate data
  pred.temp.data <- pred.temp.data[rep(1,nrow(temp.join)),]
  
  # Join both
  pred.temp.data <- cbind(pred.temp.data,temp.join)
  
  # Get LST sd and mean to rescale LST for plotting
  lst_sd <- sd(df.raw$LST_C_mean_lag)
  lst_mean <- mean(df.raw$LST_C_mean_lag)
  
  with(pred.temp.data, plot(LST_C_mean_lag*lst_sd + lst_mean, LST_C_mean_lag_squared))
  
  # Predict malaria incidence using new dataset
  mean.temp.preds <- predict(model, type = "response", newdata = pred.temp.data, se.fit = TRUE)

  
  # Make data frame with rescaled LST and lower and upper SE intervals
  marg.eff <- data.frame(LST = pred.temp.data$LST_C_mean_lag*lst_sd + lst_mean, pred = mean.temp.preds$fit, 
                         lower = mean.temp.preds$fit - mean.temp.preds$se.fit,
                         higher = mean.temp.preds$fit + mean.temp.preds$se.fit)
  
  fig2b <- ggplot(marg.eff, aes(x = LST, y = pred)) +
    geom_line(size = 1) + 
    geom_ribbon(aes(ymin = lower, ymax = higher), alpha = 0.5) +
    labs(x = "Mean Land Surface Temperature (C)", y = "Predicted malaria incidence (p.t.)") + 
    theme_minimal()
  
  fig2 <- ggarrange(fig2a, fig2b, widths = c(2, 1.7), labels = c('A','B'),font.label = list(size = 10))
  
  ggsave("figures/fig2.pdf",fig2, units="mm",width=250, height= 112, dpi=300)
  ggsave("figures/Fig2.tiff",fig2, units="in",width= 7.5, height= 3.6, dpi=300)
  
}

fig_pred <- function(model, df.norm.high.season){
  
  ### Figure 3A: Spatial predictions and residuals
  
  df.norm.high.season$pred <- predict(model, type = 'response')
  df.norm.high.season$diff <- df.norm.high.season$pred - df.norm.high.season$malaria_total_pt 
  
  high.season.ave <- aggregate(df.norm.high.season,by=list(df.norm.high.season$ID), mean, na.rm = TRUE)
  high.season.ave$ID <- as.numeric(as.character(high.season.ave$Group.1))
  
  # Load Fokontany boundary
  fkt.bound <- st_read('data/Limite_FKT/Limite_FKT_Distr_Ifanadiana.shp')
  fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))
  
  # Join averaged dataset with spatial boundary
  map.pred <- inner_join(high.season.ave, fkt.bound, by = "ID")
  
  plot.pred <- ggplot(map.pred) +
    geom_sf(aes(fill = pred, geometry = geometry))+
    scale_fill_gradient(low = 'white', high = 'royalblue4', limits = c(0,200)) + 
    labs(title = "Predicted incidence \n (per thousand)", fill = '') +
    guides(fill = guide_colourbar(barwidth = 5, barheight = 0.7)) + 
    theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=8), 
          plot.margin=margin(r=0,l=10, t=20, b=10),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank(),
          legend.position = 'bottom')

  plot.diff <- ggplot(map.pred) +
    geom_sf(aes(fill = diff, geometry = geometry))+
    labs(title = "Predicted - Observed \n (per thousand)", fill = '')  +
    scale_fill_gradient2(limits = c(-35,35)) +
    guides(fill = guide_colourbar(barwidth = 5, barheight = 0.7)) + 
    theme(plot.title = element_text(face = 'bold', hjust = 0.5,size=8), 
          plot.margin=margin(r=0,l=10, t=20, b=10),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank(),
          legend.position = 'bottom')
  
  ### Figure 3B: Within-sample predictions
  
  high.season.ave$rank.malaria <- rank(high.season.ave$malaria_total_pt)
  high.season.ave$rank.pred <- rank(high.season.ave$pred)
  
  fig3b <- ggplot(high.season.ave, aes(x = rank.malaria, y = rank.pred)) +
    geom_point(alpha = 0.7)+
    labs(x = "Malaria incidence (rank)", y = 'In-sample pred. (rank)') +
    theme_minimal() +
    theme(axis.title =element_text(size=8)) 
  
  ### Figure 3C: Out-of-sample predictions
  
  high.season.train <- df.norm.high.season[df.norm.high.season$season %in% c(1,2),]
  high.season.test <- df.norm.high.season[df.norm.high.season$season %in% c(3),]
  
  model.formula <- model$call$formula
  
  best.model.train <- glmmTMB(model.formula, data = high.season.train, ziformula = ~., family = "nbinom2",na.action="na.exclude")
  #summary(best.model.train)
  
  # Predict on test data
  high.season.test$pred <- predict(best.model.train, high.season.test, type = 'response')
  high.season.test$diff <- high.season.test$pred - high.season.test$malaria_total_pt
  
  # Calculate spatial rank
  high.season.test.av <- aggregate(high.season.test,by=list(high.season.test$ID), mean, na.rm = TRUE)
  high.season.test.av$rank.malaria <- rank(high.season.test.av$malaria_total_pt)
  high.season.test.av$rank.pred <- rank(high.season.test.av$pred)
  
  fig3c <- ggplot(high.season.test.av, aes(x = rank.malaria, y = rank.pred)) +
    geom_point(alpha = 0.7)+
    labs(x = "Malaria incidence (rank)", y = 'Out-of-sample pred. (rank)') +
    theme_minimal() +
    theme(axis.title =element_text(size=8)) 
  
  fig3b_3c <- ggarrange(fig3b, fig3c, nrow = 2, labels = c('B','C'), vjust = c(1,-0.1), font.label = list(size = 10))
  fig3 <- ggarrange(plot.pred, plot.diff, fig3b_3c, ncol = 3, labels = c('A','', ''), font.label = list(size = 10))
  fig3
  
  ggsave("figures/fig3.pdf",fig3, units="mm",width= 200, height= 150, dpi=300)
  ggsave("figures/Fig3.tiff",fig3, units="in",width= 7, height= 3.5, dpi=300)
  
  
}

fig_sem <- function(df.norm.high.season){
  
  high.season.averaged <- aggregate(df.norm.high.season,by=list(df.norm.high.season$ID), mean, na.rm = TRUE)
  
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
  
  psem_summary <- summary(psem_model_1)
  
  variable_names <- c('Malaria','Bed net','Temperature','Precipitation','Altitude','Forest Edge','Forest loss', ' Wealth score',
                      'Distance to road', 'Distance to HC', 'Residential', 'Distance to forest', ' Rice')
  
  
  # Assign width based on coefficients
  width <- 0.5 + abs(psem_summary$coefficients$Std.Estimate*5)
  
  # Assign line style based on p-values
  style <- psem_summary$coefficients$P.Value < 0.05
  style[style == TRUE] <- 'normal'
  style[style == FALSE] <- 'dashed'
  
  # Assign color based on estimate's sign
  color <- psem_summary$coefficients$Std.Estimate > 0
  color[color == TRUE] <- 'LightSeaGreen'
  color[color == FALSE] <- 'Salmon2'
  
  # Assign estimate values to arrow labels
  label_arrow <- round(psem_summary$coefficients$Std.Estimate,2)
  
  # Remove dots in variables name
  predictor_names <- gsub('[.]','',psem_summary$coefficients$Predictor)
  response_names <- gsub('[.]','',psem_summary$coefficients$Response)
  
  fig4 <- grViz("

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
  
  fig4 %>%
    export_svg %>% 
    charToRaw %>% 
    rsvg_pdf("figures/fig4.pdf")
  
}

fig_climate <- function(df.raw){
  
  malaria.per.month <- aggregate(df.raw,by=list(df.raw$ID,df.raw$month), mean, na.rm = TRUE)
  climate.per.month <- aggregate(malaria.per.month,by=list(malaria.per.month$month), mean, na.rm = TRUE)
  climate.per.month <- dplyr::select(climate.per.month, c(month,LST_C_mean_lag,Precipitation_lag))
  
  ### Panel A: Seasonal malaria incidence
  
  malaria.per.month$month <- as.factor(malaria.per.month$month)
  
  figS1A <- ggplot(malaria.per.month, aes(x = month, y = malaria_total_pt)) + 
    geom_boxplot() +
    labs(y = 'Malaria incidence (per thousand)', x = '') +
    scale_x_discrete("Month", breaks = 1:12, labels = month.abb)  +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45)) 
  
  ### Panel B: Monthly temperature and precipitation
  
  # Correct for lag
  climate.per.month$month <- climate.per.month$month - 1
  climate.per.month$month[climate.per.month$month == 0] <- 12
  
  # Define axis limits
  ylim.prim <- c(40, 372)   # in this example, precipitation
  ylim.sec <- c(21, 33) # temperature
  
  # Modifiers to plot temperature
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1] # there was a bug here
  
  figS1B <- ggplot(climate.per.month, aes(month, Precipitation_lag)) +
    geom_col(fill = 'grey') +
    geom_line(aes(y = a + LST_C_mean_lag*b), color = "red", size = 1.2) +
    scale_y_continuous("Precipitation (mm)", sec.axis = sec_axis(~ (. - a)/b, name = "Temperature (C)")) +
    scale_x_continuous("Month", breaks = 1:12, labels = month.abb)  +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45)) 
  
  figS1 <- ggarrange(figS1A, figS1B, labels = c("A","B"))

  ggsave("figures/figS1.pdf", figS1, units="mm",width=200, height= 100, dpi=300)
  
}

fig_map_variables <- function(df.averaged.raw){
  
  fkt.bound <- st_read('data/Limite_FKT/Limite_FKT_Distr_Ifanadiana.shp')
  fkt.bound <- dplyr::select(fkt.bound, c(ID, geometry))

  map.raw <- inner_join(df.averaged.raw, fkt.bound, by = "ID")
  
  # Wealth index
  plot.wealth <- ggplot(map.raw) +
    geom_sf(aes(fill = wscore.n, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Wealth index',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Distance to health center
  plot.dist <- ggplot(map.raw) +
    geom_sf(aes(fill = real.dist.csb, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Dist. to HC (km)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Bed net use
  plot.bednet <- ggplot(map.raw) +
    geom_sf(aes(fill = highseason*100, geometry = geometry), lwd = 0) +
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Bed net use (%)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Residential area
  plot.res <- ggplot(map.raw) +
    geom_sf(aes(fill = Residential*100, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Residential area (%)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Rice field area
  plot.rice <- ggplot(map.raw) +
    geom_sf(aes(fill = Rice*100, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Rice field area (%)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Distance from Residential Area to Forest
  plot.distance.forest <- ggplot(map.raw) +
    geom_sf(aes(fill = toForest_meanDist/1000, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Distance to forest (km)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Forest edge
  plot.edge.forest <- ggplot(map.raw) +
    geom_sf(aes(fill = edge_forest/1000, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Forest edge (km)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  
  # Forest loss (3y)
  plot.loss.3 <- ggplot(map.raw) +
    geom_sf(aes(fill = loss_3y*100, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Forest loss, 3y (%)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Forest loss (10y)
  plot.loss.10 <- ggplot(map.raw) +
    geom_sf(aes(fill = loss_10y*100, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Forest loss, 10y (%)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Minimum LST
  plot.temp.min <- ggplot(map.raw) +
    geom_sf(aes(fill = LST_C_min_lag, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Min. LST (C)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Maximum LST
  plot.temp.max <- ggplot(map.raw) +
    geom_sf(aes(fill = LST_C_max_lag, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Max. LST (C)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Mean LST
  plot.temp.mean <- ggplot(map.raw) +
    geom_sf(aes(fill = LST_C_mean_lag, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Mean LST (C)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Suitability index
  plot.temp.mean.squared <- ggplot(map.raw) +
    geom_sf(aes(fill = LST_C_mean_lag_squared, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Suitability index',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  
  # Precipitation
  plot.prec<- ggplot(map.raw) +
    geom_sf(aes(fill = Precipitation, geometry = geometry), lwd = 0)+
    scale_fill_gradient(low = 'gray96', high = 'royalblue4') + 
    labs(title = 'Monthly precipitation (mm)',fill = '') +
    theme(legend.title=element_text(size=7),legend.text=element_text(size=8),
          plot.title = element_text(size = 9, hjust = 0.5), legend.position = 'right',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank()) + 
    guides(fill = guide_colourbar(barwidth = 0.7, barheight = 5)) 
  

  figs2<- ggarrange(plot.wealth, plot.dist, plot.bednet, plot.res, plot.rice,  plot.distance.forest,
                            plot.edge.forest, plot.loss.10, plot.loss.3, 
                            plot.temp.min, plot.temp.max, plot.temp.mean, plot.temp.mean.squared, plot.prec)
  
  ggsave("figures/figs2.pdf", figs2, units="mm",width=200, height= 170, dpi=300)
  
  
}

fig_malaria_seasonal <- function(model, df.norm.high.season){
  
  df.norm.high.season$pred <- predict(model, type = 'response')
  df.norm.high.season$diff <- df.norm.high.season$malaria_total_pt - df.norm.high.season$pred
  
  df.norm.high.season$time.labels <- paste(month.abb[df.norm.high.season$month],df.norm.high.season$year)
  df.norm.high.season$time.labels <- factor(df.norm.high.season$time.labels, levels = c("Apr 2014", "May 2014","Jun 2014", "Jul 2014", "Aug 2014", "Sep 2014", "Oct 2014", "Nov 2014", "Dec 2014", "Jan 2015", "Feb 2015",
                                                                          "Mar 2015", "Apr 2015", "May 2015", "Jun 2015", "Jul 2015", "Aug 2015", "Sep 2015", "Oct 2015", "Nov 2015", "Dec 2015",
                                                                          "Jan 2016", "Feb 2016", "Mar 2016", "Apr 2016", "May 2016", "Jun 2016", "Jul 2016", "Aug 2016", "Sep 2016", "Oct 2016",
                                                                          "Nov 2016", "Dec 2016", "Jan 2017", "Feb 2017", "Mar 2017", "Apr 2017", "May 2017", "Jun 2017", "Jul 2017", "Aug 2017",
                                                                          "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017"))
  
  
  high.season.pred.time <- aggregate(df.norm.high.season,by=list(df.norm.high.season$time), mean, na.rm = TRUE)
  
  small.pred <- dplyr::select(df.norm.high.season, c('ID','month','year','time.labels','malaria_total_pt','pred'))
  long.pred <- gather(small.pred, group, value, c('malaria_total_pt','pred'), factor_key=TRUE)
  #long.pred$time <- as.factor(long.pred$time)
  
  figs3 <- ggplot(long.pred, aes(x=time.labels, y=value, fill=group)) + 
    geom_boxplot(outlier.shape = NA) +
    labs(y = 'Malaria incidence (per thousand)', x = '', fill = '') +
    scale_fill_discrete(labels = c("Observed", "Predicted")) +
    scale_y_continuous(limits = c(0,400)) +
    theme_minimal()  +
    theme(legend.position="bottom", axis.text.x=element_text(angle = 45, hjust = 1), plot.margin=margin(r=10,l=30, t=20, b=10), axis.title=element_text(size=9))
    
  ggsave("figures/figs3.pdf", figs3, units="mm",width=200, height= 100, dpi=300)
  
}

#### 2. RUN ANALYSIS ####

# Load all datasets

malaria_raw <- read.csv('output/malaria_v3.csv')
averaged.raw <- read.csv('output/averaged_raw.csv')
malaria.norm <- read.csv('output/malaria_norm.csv')
malaria.norm.high.season <- read.csv('output/malaria_norm_high_season.csv')
high.season.averaged <- read.csv('output/high_season_averaged.csv')


# Note: I retrain models here as saved models cannot always be read if glmmTMB is updated. 

# Final processing for the covariance structures
malaria.norm.high.season$time_factor <- numFactor(malaria.norm.high.season$time)

malaria.norm.high.season$pos <- numFactor(scale(malaria.norm.high.season$long), scale(malaria.norm.high.season$lat))
malaria.norm.high.season$group <- 'All'

# Train model with month and Fokontany as random effects
re.model <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                       Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                       LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + 
                       (1|ID) + (1|month), 
                     data = malaria.norm.high.season, 
                     family = "nbinom2", ziformula = ~., na.action = "na.exclude")
summary(re.model)

# Train same model without random effects to compare
fixed.model <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                      Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                      LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared, 
                    data = malaria.norm.high.season, 
                    family = "nbinom2", ziformula = ~., na.action = "na.exclude")


# Train model that explicitly accounts for temporal and spatial auto-correlation 

sp.temp.model <- glmmTMB(formula = malaria_total_pt ~ highseason + Residential + 
                      Rice + real.dist.csb + wscore.n + loss_3y + Precipitation_lag + 
                      LST_C_mean_lag + LST_C_min_lag + LST_C_mean_lag_squared + 
                        ou(time_factor-1 |Commune) + mat(pos + 0|group), 
                    data = malaria.norm.high.season, 
                    family = "nbinom2", ziformula = ~., na.action = "na.exclude")
summary(sp.temp.model)

# Note that this model cannot easily be used for forward predictions, which is why I use re.model for predictions

### 2.1 Analysis of malaria incidence and other variables

out_variable_analysis <- analysis_variables(malaria_raw, averaged.raw)

malaria_metrics <- out_variable_analysis$malaria
malaria_metrics

monthly_data <- out_variable_analysis$climate
monthly_data

### 2.2 Check for spatial and temporal autocorrelation in both models

# Compare temporal autocorrelation between RE vs OU/Matern model
re.acf <- calculateACF(re.model, malaria.norm.high.season)
re.acf$test$p.value

sp.temp.acf <- calculateACF(sp.temp.model, malaria.norm.high.season)
sp.temp.acf$test$p.value

# Compare spatial autocorrelation between RE vs OU/Matern model
re.moran <- calculateMoran(re.model, malaria.norm.high.season)
re.moran$p.value
re.moran$observed

sp.temp.moran <- calculateMoran(sp.temp.model, malaria.norm.high.season)
sp.temp.moran$p.value
sp.temp.moran$observed

### 2.3 Analysis of model predictions

# Stats used in the text
in_sample_pred(re.model, malaria.norm.high.season)
out_sample_pred(re.model, malaria.norm.high.season)

in_sample_pred(fixed.model, malaria.norm.high.season)
out_sample_pred(fixed.model, malaria.norm.high.season)

#### 3. MAKE TABLES ####

### Table 1: Model variables 

# See text

### Table S1: Model coefficients for model with OU and Matern

table_model_coefficients(sp.temp.model)

### Table S2: Model selection steps

# See text

### Table S3: Model average 

table_model_avg()

### Table S4: Model coefficients for model with random effects

table_model_coefficients(re.model)

###

#### 4. MAKE FIGURES ####

### Figure 1: Malaria and land-use maps

fig_map_malaria(averaged.raw)

### Figure 2: Model coefficients and marginal effect of LST

fig_coefficients(sp.temp.model, malaria.norm.high.season, malaria_raw)

### Figure 3: Model predictions

fig_pred(re.model, malaria.norm.high.season)

### Figure 4: SEM 

fig_sem(malaria.norm.high.season)

### Figure S1: Seasonal malaria patterns and climate

fig_climate(malaria_raw)

### Figure S2: Maps for all variables

fig_map_variables(averaged.raw)

### Figure S3: Seasonal malaria incidence and predictions

fig_malaria_seasonal(re.model, malaria.norm.high.season)





