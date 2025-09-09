## MESS
# generate multivariate environmental similarity surface using
# dismo, which works only with raster :eyeroll: 
# per Elith et al. 2010
# https://doi.org/10.1111/j.2041-210X.2010.00036.x


library(terra)
library(geodata)
library(tidyverse)


# get administrative area for a country
somalia_shp <- gadm(
   country = "SOM",
   level = 0,
   path = tempdir()
)

# have a look at it
plot(somalia_shp)

# download climactic data for this country
bioclim_somalia <- worldclim_country(
  country = "SOM",
  var = "bio",
  res = 0.5,
  path = "data/downloads"
)

# subset to just
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO12 = Annual Precipitation
covs_somalia <- bioclim_somalia[[c(5, 6, 12)]] |>
  # and mask to country shapefile
  mask(somalia_shp)

names(covs_somalia) <- c("tmax_warm", "tmin_cool", "precip")

# have a look
plot(covs_somalia)


# make up some points that fall within the boundary of your
# shapefile
somalia_pts <- tibble(
  x = c(43, 45, 44, 44, 45, 50),
  y =c(2, 3, 4, 2.1, 4, 10)
)

# check that the points fall inside and edit above until they do
plot(somalia_shp)
points(somalia_pts)


# see points against each of the layers
par(mfrow = c(2,2))
plot(covs_somalia[[1]])
points(somalia_pts)

plot(covs_somalia[[2]])
points(somalia_pts)

plot(covs_somalia[[3]])
points(somalia_pts)

par(mfrow = c(1,1))

#####
# turn covariates into raster package format for mess
library(dismo)
library(raster)

# extract covariate values at our points
coord_covs <- extract(
  covs_somalia,
  somalia_pts
)

# convert layers to raster package for mess
covraster <- raster::brick(covs_somalia)

# calculate multivariate environmental similarity surface (MESS)
# for our points against whole raster
mess_somalia <- mess(
  x = covraster,
  v = coord_covs |>
    select(-ID) |>
    as.data.frame()
) |>
  rast()

# plot the result
plot(mess_somalia)
points(somalia_pts)
# this result is in a unitless format but the closer to zero
# the value is, the more similar

par(mfrow = c(2,2))
plot(covs_somalia[[1]])
points(somalia_pts)

plot(covs_somalia[[2]])
points(somalia_pts)

plot(covs_somalia[[3]])
points(somalia_pts)

plot(mess_somalia)
points(somalia_pts)

par(mfrow = c(1,1))


# plot with palette that diverges around zero
plot_mess_local <- ggplot() +
  geom_spatraster(
    data = mess_somalia
  ) +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = 1,
    limit = max(abs(values(mess_somalia)), na.rm = TRUE) * c(-1, 1)
  ) +
  theme_void() +
  labs(fill = "Multivariate\nEnvironmental\nSimilarity")

plot_mess_local





# make mask of this, such that anything < 0 is NA,
# i.e. dissimilar, and >= 0 is 1, i.e., similar.
mess_mask <- mess_somalia

mvals <- values(mess_somalia)

mess_mask[which(mvals < 0)] <- NA
mess_mask[which(mvals >= 0)] <- 1

mess_mask <- mask(mess_mask, somalia_shp)

# look at area of country we have represented
plot(mess_mask)
plot(somalia_shp, add = TRUE)
points(somalia_pts)

# check the covariates we have represented
masked_covs <- mask(covs_somalia, mess_mask)
plot(masked_covs)


