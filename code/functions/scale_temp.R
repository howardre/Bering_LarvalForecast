scale_temp <- function(gam, df){
  values <- visreg::visreg(gam[[2]], 
                           data = df,
                           type = "contrast", 
                           plot = FALSE)
  smooths <- plyr::ldply(values, function(part)
    data.frame(variable = part$meta$x,
               x = part$fit[[part$meta$x]],
               smooth = part$fit$visregFit,
               lower = part$fit$visregLwr,
               upper = part$fit$visregUpr))
  temp_smooth <- smooths %>%
    filter(variable == "roms_temperature" &
             lower > 0)
  return(temp_smooth)
}

scale_temp1 <- function(gam, df){
  values <- visreg::visreg(gam[[2]], 
                           data = df,
                           type = "contrast", 
                           plot = FALSE)
  smooths <- plyr::ldply(values, function(part)
    data.frame(variable = part$meta$x,
               x = part$fit[[part$meta$x]],
               smooth = part$fit$visregFit,
               lower = part$fit$visregLwr,
               upper = part$fit$visregUpr))
  temp_smooth <- smooths %>%
    filter(variable == "roms_temperature" &
             smooth > 0)
  return(temp_smooth)
} # only use this if there are no values within range
