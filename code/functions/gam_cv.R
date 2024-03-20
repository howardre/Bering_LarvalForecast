gam_cv <- function(data){
  train <- data %>% filter(year < 2011)
  test <- data %>% filter(year >= 2011)
  gam_small <- gam(larvalcatchper10m2 ~ s(year, bs = 're') +
                     te(lon, lat, bs = "tp", m = 1) + 
                     s(doy, k = 9, bs = "tp", m = 1),
                   data = train,
                   family = tw(link = 'log'),
                   method = 'REML')
  gam_base <- gam(larvalcatchper10m2 ~ s(year, bs = 're') + 
                    s(doy, k = 9, bs = "tp", m = 1) +
                    te(lon, lat, bs = "tp", m = 1) +
                    s(roms_temperature, k = 5, bs = "tp", m = 1) +
                    s(roms_salinity, k = 5, bs = "tp", m = 1),
                  data = train,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(larvalcatchper10m2 ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "tp", m = 1) +
                       te(lon, lat, bs = "tp", m = 1) +
                       s(roms_temperature, k = 5, bs = "tp", m = 1) +
                       s(roms_salinity, k = 5, bs = "tp", m = 1) +
                       s(doy, by = mean_temp, bs = "tp", m = 1),
                     data = train,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog <- gam(larvalcatchper10m2 ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "tp", m = 1) +
                       te(lon, lat, bs = "tp", m = 1) +
                       s(roms_temperature, k = 5, bs = "tp", m = 1) +
                       s(roms_salinity, k = 5, bs = "tp", m = 1) +
                       s(lon, lat, by = mean_temp, bs = "tp", m = 1),
                     data = train,
                     family = tw(link = 'log'),
                     method = 'REML')
  
  # Get correlation coefficient
  test$small_pred <- predict(gam_small,
                             newdata = test,
                             type = "response")
  small_sp <- cor.test(test$larvalcatchper10m2, 
                       test$small_pred, 
                       method = 'spearman',
                       exact = FALSE)
  test$base_pred <- predict(gam_base,
                            newdata = test,
                            type =  "response")
  base_sp <- cor.test(test$larvalcatchper10m2, 
                      test$base_pred, 
                      method = 'spearman',
                      exact = FALSE)
  test$vc_phen_pred <- predict(gam_vc_phen,
                               newdata = test,
                               type = "response")
  vc_phen_sp <- cor.test(test$larvalcatchper10m2,
                         test$vc_phen_pred,
                         method = 'spearman',
                         exact = FALSE)
  test$vc_geog_pred <- predict(gam_vc_geog,
                               newdata = test,
                               type = "response")
  vc_geog_sp <- cor.test(test$larvalcatchper10m2,
                         test$vc_geog_pred,
                         method = 'spearman',
                         exact = FALSE)

  gam_list <- list(small_sp, base_sp, vc_phen_sp, vc_geog_sp)
}
