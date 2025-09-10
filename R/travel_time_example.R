library(traveltime)
library(tidyverse)
library(terra)
library(geodata)






# get administrative area for a country
somalia_shp <- gadm(
  country = "SOM",
  level = 0,
  path = tempdir()
)

# have a look at it
plot(somalia_shp)

# make up some points that fall within the boundary of your
# shapefile
somalia_pts <- tibble(
  x = c(43.7696, 45.3182, 44.0536, 44.7721, 45.3496, 50),
  y =c(2, 3, 4, 2.1, 4, 10)
)

# check that the points fall inside and edit above until they do
plot(somalia_shp)
points(somalia_pts)

# download friction surface for our shapefile
somalia_friction <- get_friction_surface(
  surface = "motor2020",
  filename = tempfile(fileext = ".tif"),
  extent = somalia_shp
)

plot(somalia_friction)

somalia_friction_masked <- somalia_friction |>
  mask(mask = somalia_shp)

mask(
  x = somalia_friction,
  mask = somalia_shp
)


plot(somalia_friction_masked)

