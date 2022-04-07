gam_selection <- function(data){
  gam_base <- gam(log(larvalcatchper10m2 + 1) ~ s(year, bs = 're') + 
                    s(doy, k = 8) +
                    s(lon, lat) +
                    s(roms_temperature, k = 6) +
                    s(roms_salinity, k = 6),
                  data = data[data$larvalcatchper10m2 > 0, ])
  gam_vc_phen <- gam(log(larvalcatchper10m2 + 1) ~ s(year, bs = 're') +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(doy, by = mean_temp, k = 6),
                     data = data[data$larvalcatchper10m2 > 0, ])
  gam_vc_geog <- gam(log(larvalcatchper10m2 + 1) ~ s(year, bs = 're') +
                       s(doy, k = 8) +
                       s(lon, lat) +
                       s(roms_temperature, k = 6) +
                       s(roms_salinity, k = 6) +
                       s(lon, lat, by = mean_temp, k = 6),
                     data = data[data$larvalcatchper10m2 > 0, ])
  gam_list <- list(gam_base, gam_vc_phen, gam_vc_geog)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}
