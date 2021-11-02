gam_selection <- function(data, roms_temps){
  gam_base <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                    s(doy, k = 8) +
                    s(lon, lat) +
                    s(roms_temperature, k = 6) +
                    s(roms_salinity, k = 6),
                  data = data,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(doy, by = mean_temp, k = 6),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_interaction <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                           s(doy, k = 8) +
                           s(lon, lat) +
                           s(roms_temperature, roms_salinity),
                         data = data,
                         family = tw(link = 'log'),
                         method = 'REML')
  gam_interaction2 <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                            s(doy, k = 8) +
                            s(lon, lat) +
                            s(roms_temperature, roms_salinity) +
                            s(roms_temperature, k = 6) +
                            s(roms_salinity, k = 6),
                          data = data,
                          family = tw(link = 'log'),
                          method = 'REML')
  gam_vc_phen2 <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, roms_salinity) +
                        s(roms_temperature) +
                        s(roms_salinity, k = 4) +
                        s(doy, by = mean_temp, k = 6),
                      data = data,
                      family = tw(link = 'log'),
                      method = 'REML')
  gam_vc_geog <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(lon, lat, by = mean_temp, k = 6),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog2 <- gam(larvalcatchper10m2 + 1 ~ factor(year) +
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, roms_salinity) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(lon, lat, by = mean_temp, k = 6),
                      data = data,
                      family = tw(link = 'log'),
                      method = 'REML')
  gam_list <- list(gam_base, gam_interaction, gam_interaction2, gam_vc_phen, 
                   gam_vc_phen2, gam_vc_geog, gam_vc_geog2)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}