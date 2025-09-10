# HGAM example


library(terra)
source("R/functions.R")

# should run without part 1 if data/grids/ contains the below rasters.
# as an alternative to running script preparing_data_1.R, grids can be
# downloaded into data/grids from:
# https://doi.org/10.26188/21630146.v1 
# you can do this manually or work out how to in R!


# load in real data
kenya_mask <- terra::rast("data/grids/kenya_mask.tif")
bc_kenya <- terra::rast("data/grids/bc_kenya.tif")
rescale_travel <- terra::rast("data/grids/rescale_travel.tif")


#specify covariates
# BIO4 = Temperature Seasonality (standard deviation ×100)
tseas <- bc_kenya[[4]] 
# BIO5 = Max Temperature of Warmest Month
tmax   <- bc_kenya[[5]]
# BIO7 = Temperature Annual Range (BIO5-BIO6)
trange <- bc_kenya[[7]]

covs <- c(tseas, tmax, trange)
names(covs) <- c("tseas", "tmax", "trange")

# bias layer
bias <- rescale_travel ^ 2
names(bias) <- "bias"

# terra::writeRaster(
#   x = bias,
#   filename = "data/grids/bias.tif"
# )

# generate the fake 'true' relative abundance raster

# generate an unscaled relative abundance
# rel_abund_unscaled <- exp(-1 + covs$tmax * 0.1)
# 
# plot(rel_abund_unscaled)

# create your own here:

ra_complex <- exp(
  1 +
    covs$tseas * 0.01 +
    covs$tmax * -0.03 +
    covs$trange * -0.01 + 
    covs$trange ^ 2 * 0.02
)

ra_sp_1 <- ra_complex + exp(
  covs$tseas^2 * 0.00005 +
    covs$tmax * -0.01 +
    covs$trange ^ 2 * 0.05
)


ra_sp_2 <- ra_complex + exp(
  covs$tseas * 0.005 +
    covs$trange * 0.01 +
    covs$tmax ^ 2 * 0.001
)

plot(c(ra_complex, ra_sp_1, ra_sp_2))



log(abund[i]) ~ group_intercept + group_coef_temp*temp + species_coef_temp[i]*temp

species_coef_temp[3] <- 0

species_coef_temp[1:2] <- rnorm(2, 0, 1)