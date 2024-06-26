# Libraries
library(ggOceanMaps)
library(here)

# Create palettes and data frames
bathy_palette <- colorRampPalette(c("lightsteelblue1", "lightsteelblue4"))(9)
BS_bathy <- data.frame(lon = c(-176.6, -158.5),
                       lat = c(51.5, 62.5))
text_labels <- data.frame(place = c("Pribilof Islands", "St. Matthew Island", 
                                    "Alaska Peninsula", "Unimak Island", 
                                    "Unimak Pass", "Unalaska Island"),
                          lat = c(56.5, 60, 57, 54.2, 54.6, 53.8),
                          lon = c(-172.6, -172, -162.2, -162.4, -167, -169.5))
text_sf <- sf::st_as_sf(text_labels, 
                        coords = c("lon", "lat"),
                        crs = sf::st_crs(4326))

alaska_label <- data.frame(place = "Alaska",
                           lat = 60.5,
                           lon = -158)
alaska_sf <- sf::st_as_sf(alaska_label, 
                          coords = c("lon", "lat"),
                          crs = sf::st_crs(4326))

# Make map
bering_map <- basemap(data = BS_bathy,
                      bathymetry = TRUE,
                      rotate = TRUE,
                      legends = FALSE,
                      land.col = "wheat4",
                      grid.col = NA,
                      lon.interval = 5,
                      lat.interval = 2) +
  geom_polygon(data = transform_coord(BS_bathy),
               aes(x = lon, y = lat),
               fill = NA) +
  scale_fill_manual(values = bathy_palette) +
  labs(x = "Longitude",
       y = "Latitude")  +
  ggspatial::annotation_north_arrow(location = "tr",
                                    which_north = "true",
                                    style = ggspatial::north_arrow_nautical(text_family = "serif"),
                                    height = unit(3, "cm"),
                                    width = unit(3, "cm")) +
  ggspatial::annotation_scale(location = "bl",
                              text_family = "serif",
                              style = "ticks",
                              height = unit(0.7, "cm"),
                              text_cex = 1.6,
                              line_width = 1.8) +
  geom_sf_text(data = text_sf,
               aes(label = place),
               size = 7,
               family = "serif") +
  geom_sf_text(data = alaska_sf,
               aes(label = place),
               size = 9,
               family = "serif",
               fontface = "bold") +
  theme(axis.text = element_text(family = "serif", size = 20),
        axis.title = element_text(family = "serif", size = 23),
        strip.text = element_text(family = "serif", size = 21))

bering_map

# Save map
bering_map 
dev.copy(jpeg,
         here('manuscripts/Figures',
              'Figure1.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Presentation map
pk_larvae <- readRDS(here('data', 'pk_larvae.rds'))

basemap(data = BS_bathy,
        bathymetry = TRUE,
        rotate = TRUE,
        legends = FALSE,
        land.col = "wheat4",
        grid.col = NA,
        lon.interval = 5,
        lat.interval = 2) +
  geom_polygon(data = transform_coord(BS_bathy),
               aes(x = lon, y = lat),
               fill = NA) +
  scale_fill_manual(values = bathy_palette) +
  labs(x = "Longitude",
       y = "Latitude")  +
  ggspatial::annotation_north_arrow(location = "tr",
                                    which_north = "true",
                                    style = ggspatial::north_arrow_nautical(text_family = "sans"),
                                    height = unit(3, "cm"),
                                    width = unit(3, "cm")) +
  ggspatial::annotation_scale(location = "bl",
                              text_family = "sans",
                              style = "ticks",
                              height = unit(0.7, "cm"),
                              text_cex = 1.6,
                              line_width = 1.8) +
  ggspatial::geom_spatial_point(data = pk_larvae,
                                aes(x = lon, y = lat), 
                                color = "firebrick4",
                                size = 2.5) +
  theme(axis.text = element_text(family = "sans", size = 20),
        axis.title = element_text(family = "sans", size = 23),
        strip.text = element_text(family = "sans", size = 21))
dev.copy(jpeg,
         here('results',
              'pklarvae_map.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()