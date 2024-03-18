gam_selection <- function(data){
  gam_base <- gam(larvalcatchper10m2 ~ s(year, bs = 're') + 
                    s(doy, k = 9, bs = "tp", m = 1) +
                    te(lon, lat, bs = "tp", m = 1) +
                    s(roms_temperature, k = 5, bs = "tp", m = 1) +
                    s(roms_salinity, k = 5, bs = "tp", m = 1),
                  data = data,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(larvalcatchper10m2 ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "tp", m = 1) +
                       te(lon, lat, bs = "tp", m = 1) +
                       s(roms_temperature, k = 5, bs = "tp", m = 1) +
                       s(roms_salinity, k = 5, bs = "tp", m = 1) +
                       s(doy, by = mean_temp, bs = "tp", m = 1),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog <- gam(larvalcatchper10m2 ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "tp", m = 1) +
                       te(lon, lat, bs = "tp", m = 1) +
                       s(roms_temperature, k = 5, bs = "tp", m = 1) +
                       s(roms_salinity, k = 5, bs = "tp", m = 1) +
                       s(lon, lat, by = mean_temp, bs = "tp", m = 1),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_list <- list(gam_base, gam_vc_phen, gam_vc_geog)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}