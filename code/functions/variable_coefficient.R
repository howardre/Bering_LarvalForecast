variable_coefficient <- function(gam, data){
  preds <- predict(gam[[2]], type = 'terms', se.fit = T)
  pred_slope <- preds[[1]][, 6] / data$mean_temp
  pred_slope_se <- 1.96 * preds[[2]][, 6]
  pred_slope_up <- (preds[[1]][, 6] + pred_slope_se) / data$mean_temp
  pred_slope_low <- (preds[[1]][, 6] - pred_slope_se) / data$mean_temp
  sign_slope_pos <- (1:length(pred_slope))[pred_slope_low > 0]
  sign_slope_neg <- (1:length(pred_slope))[pred_slope_up < 0]
  return(list(sign_slope_neg, sign_slope_pos, pred_slope))
}
