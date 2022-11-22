gam_selection <- function(data){
  gam_base <- gam(catch ~ s(year, bs = 're') + 
                    s(doy, k = 9, bs = "gp") +
                    s(lon, lat, bs = "gp") +
                    s(roms_temperature, k = 9, bs = "gp") +
                    s(roms_salinity, k = 9, bs = "gp"),
                  data = data,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(catch ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "gp") +
                       s(lon, lat, bs = "gp") +
                       s(roms_temperature, k = 9, bs = "gp") +
                       s(roms_salinity, k = 9, bs = "gp") +
                       s(doy, by = mean_temp, bs = "gp"),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog <- gam(catch ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "gp") +
                       s(lon, lat, bs = "gp") +
                       s(roms_temperature, k = 9, bs = "gp") +
                       s(roms_salinity, k = 9, bs = "gp") +
                       s(lon, lat, by = mean_temp, bs = "gp"),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_list <- list(gam_base, gam_vc_phen, gam_vc_geog)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}
