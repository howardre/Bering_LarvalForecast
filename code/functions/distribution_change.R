distribution_change <- function(data, grids){
  my_color = colorRampPalette(c("#1C0D51", "#4C408E", "#7E77B0",
                                "#AFABCB", "#DAD9E5", "#F9F9F9",
                                "#FFDAB7", "#FFB377","#E18811",
                                "#AC6000", "#743700"))
  color_levels = 100
  max_absolute_value = max(abs(c(min(grids[[1]]$diff, na.rm = T), max(grids[[1]]$diff, na.rm = T)))) 
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(grids[[1]]$diff, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  nlat = 70
  nlon = 90
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  image(lond,
        latd,
        t(matrix(grids[[1]]$diff,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include], 
        breaks = color_sequence[breaks_to_include],
        ylab = "Longitude",
        xlab = "Latitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        main = '',
        cex.lab = 6.5,
        cex.axis = 6.5)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
}
