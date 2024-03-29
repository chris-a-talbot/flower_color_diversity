# Load the necessary libraries
library(ggplot2)
library(maps)

# Define the boundaries of the map
lon_min <- -92.5
lon_max <- -67.5
lat_min <- 38
lat_max <- 48

# Create the base map
world_map <- map_data("world")

# Plot the map with ggplot2
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
  coord_fixed(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), ratio = 1) +
  labs(title = "Zoomed-in World Map", x = "Longitude", y = "Latitude") +
  theme_minimal()
