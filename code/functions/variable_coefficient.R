variable_coefficient <- function(gam, data){
  preds <- predict(gam[[2]], type = 'terms', se.fit = T)
  data$pred_slope <- preds[[1]][, 6] / data$mean_temp
  data$pred_slope_se <- 1.96 * preds[[2]][, 6]
  data$pred_slope_up <- (preds[[1]][, 6] + data$pred_slope_se) / data$mean_temp
  data$pred_slope_low <- (preds[[1]][, 6] - data$pred_slope_se) / data$mean_temp
  data <- data %>% mutate(sign_slope_pos = ifelse(pred_slope_low > 0, pred_slope, 0))
  data <- data %>% mutate(sign_slope_neg = ifelse(pred_slope_up < 0, pred_slope, 0))
  return(data)
}
