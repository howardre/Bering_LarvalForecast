gam_selection <- function(data, roms_temps){
  gam_base <- gam(catch ~ factor(year) + 
                    s(doy, k = 8) +
                    s(lon, lat) +
                    s(roms_temperature, k = 6) +
                    s(roms_salinity, k = 6),
                  data = data,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(catch ~ factor(year) +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(doy, by = mean_temp, k = 6),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog <- gam(catch ~ factor(year) +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(lon, lat, by = mean_temp, k = 6),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_list <- list(gam_base, gam_vc_phen, gam_vc_geog)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}
