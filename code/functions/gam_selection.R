gam_selection <- function(data){
  gam_base <- gam(catch ~ s(year, bs = 're') + 
                    s(doy, k = 9, bs = "gp") +
                    te(lon, lat, bs = "gp", m = mp, d = 2) +
                    s(roms_temperature, k = 5, bs = "gp", m = c(-2, r, 0.5)) +
                    s(roms_salinity, bs = "gp", m = c(-2, r, 0.5)),
                  data = data,
                  family = tw(link = 'log'),
                  method = 'REML')
  gam_vc_phen <- gam(catch ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "gp", m = c(-2, r, 0.5)) +
                       te(lon, lat, bs = "gp", m = mp, d = 2) +
                       s(roms_temperature, k = 5, bs = "gp", m = c(-2, r, 0.5)) +
                       s(roms_salinity, k = 5, bs = "gp", m = c(-2, r, 0.5)) +
                       s(doy, by = mean_temp, bs = "gp", m = c(-2, r, 0.5)),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_vc_geog <- gam(catch ~ s(year, bs = 're') +
                       s(doy, k = 9, bs = "gp", m = c(-2, r, 0.5)) +
                       te(lon, lat, bs = "gp", m = mp, d = 2) +
                       s(roms_temperature, k = 5, bs = "gp", m = c(-2, r, 0.5)) +
                       s(roms_salinity, k = 5, bs = "gp", m = c(-2, r, 0.5)) +
                       te(lon, lat, by = mean_temp, bs = "gp", m = mp, d = 2),
                     data = data,
                     family = tw(link = 'log'),
                     method = 'REML')
  gam_list <- list(gam_base, gam_vc_phen, gam_vc_geog)
  best_gam <- gam_list[[which.min(sapply(1:length(gam_list),
                                         function(x) min(gam_list[[x]]$aic)))]]
  return_list <- list(gam_list, best_gam)
}
