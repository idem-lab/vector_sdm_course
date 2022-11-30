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

terra::writeRaster(
  x = bias,
  filename = "data/grids/bias.tif"
)

# generate the fake 'true' relative abundance raster

# generate an unscaled relative abundance
# rel_abund_unscaled <- exp(-1 + covs$tmax * 0.1)
# 
# plot(rel_abund_unscaled)

# create your own here:

rel_abund_unscaled <- exp(1 +
                            covs$tseas * -0.01 +
                            covs$tmax * -0.03 +
                            covs$trange * -0.01 + 
                            covs$trange ^ 2 * 0.02)

# rescale the relative abundance, from 0 to 1
rel_abund <- rescale_abundance(rel_abund_unscaled)

names(rel_abund) <- "relative_abundance"

plot(rel_abund)

terra::writeRaster(
  x = rel_abund,
  filename = "data/grids/rel_abund.tif"
)

# sample abundance data at a random set of locations
n_samples <- 100

# random locations all over the country - unweighted sampling
sample_locations_random <- random_locations(kenya_mask,
                                            n_samples,
                                            weighted = FALSE)

plot(rel_abund)
points(sample_locations_random)

catches_random <- sim_catches(sample_locations = sample_locations_random,
                              relative_abundance = rel_abund,
                              max_average_catch_size = 100)

plot(rel_abund)
points(catches_random, pch = 21, bg = catches_random$presence)

# random locations, biased as per the bias layer - e.g. convenience samples
sample_locations_bias_weighted <- random_locations(bias,
                                                   n_samples)
plot(bias)
points(sample_locations_bias_weighted, pch = 16)
catches_bias_weighted <- sim_catches(sample_locations = sample_locations_bias_weighted,
                                     relative_abundance = rel_abund,
                                     max_average_catch_size = 100)

plot(rel_abund)
points(catches_bias_weighted, pch = 21, bg = catches_bias_weighted$presence)

 

# random locations, biased as per the relative abundance layer - e.g. targeted
# to areas of high abundance (where malaria interventions happen?)
sample_locations_abundance_weighted <- random_locations(rel_abund ^ (1/3),
                                                        n_samples)

catches_abundance_weighted <- sim_catches(sample_locations = sample_locations_abundance_weighted,
                                          relative_abundance = rel_abund,
                                          max_average_catch_size = 100)

plot(rel_abund)
points(catches_abundance_weighted,
       pch = 21,
       bg = catches_abundance_weighted$presence)


# simulating biased occurrence data - there are two ways:

# 1. simulate biased sampling locations, and then sample the presence at those
# locations, and only keep the ones in which they are present

# filter out biased ones to only the 1s, and store as coordinates of occurrence records
sample_locations_bias_weighted <- random_locations(bias,
                                                   n_samples)
catches_bias_weighted <- sim_catches(sample_locations = sample_locations_bias_weighted,
                                     relative_abundance = rel_abund,
                                     max_average_catch_size = 100)
occurrence_coords <- crds(catches_bias_weighted[catches_bias_weighted$presence == 1])

plot(kenya_mask)
points(occurrence_coords, pch = 16)

# this is great, but when simulating, it's difficult to get a realistic number
# of occurrence records for model fitting! You can try changing the number of
# sampling locations (n_samples in the code above) until it looks good:

sample_locations_bias_weighted <- random_locations(bias,
                                                   200)  # <- change this number
catches_bias_weighted <- sim_catches(sample_locations = sample_locations_bias_weighted,
                                     relative_abundance = rel_abund,
                                     max_average_catch_size = 100)
occurrence_coords <- crds(catches_bias_weighted[catches_bias_weighted$presence == 1])

# number of records
nrow(occurrence_coords)
plot(kenya_mask)
points(occurrence_coords, pch = 16)


# 2. alternatively, you can simulate the distribution of occurrence points by
# simulating locations, biased by the *product* of the bias and probability of
# detection in each cell. The probability of detection is calculated using a
# formula that assumes the occurrence comes from a catch (like the simulations
# above), and that the number of mosquitoes caught is a Poisson sample, given
# the average number. It is the probability of observing one *or more*
# mosquitoes in the catch. 

prob_present <- probability_of_presence(rel_abund, max_average_catch_size = 100)
names(prob_present) <- "prob_present"
plot(prob_present)

reported_occurrence_rate <- bias * prob_present
names(reported_occurrence_rate) <- "rep_occ_rate"

plot(reported_occurrence_rate)
n_occurrences <- 100
sample_locations_bias_weighted <- random_locations(reported_occurrence_rate,
                                                   n_occurrences)

# plot the reported occurrence rates
plot(reported_occurrence_rate)
points(sample_locations_bias_weighted, pch = 16)

plot(
  c(
    rel_abund,
    prob_present,
    reported_occurrence_rate,
    bias
  )
)

terra::writeRaster(
  x = prob_present,
  filename = "data/grids/prob_present.tif"
)

terra::writeRaster(
  x = reported_occurrence_rate,
  filename = "data/grids/reported_occurrence_rate.tif"
)

dir.create("data/tabular")
# save occurrence data

# presence_only
write.csv(
  x = occurrence_coords,
  file = "data/tabular/presence_only_data.csv",
  row.names = FALSE
)

# format presence absence data

# presence-absence data with sampling locations randomly selected
pa_random_data <- as_tibble(catches_random) %>%
  dplyr::select(count, presence) %>%
  bind_cols(as_tibble(crds(catches_random)))

# # non-tidyverse way
# pa_random_occurrence <- as.data.frame(catches_random)
# pa_random_occurrence <- pa_random_occurrence[, c("count", "presence")]
# pa_random_coords <- as.data.frame(crds(catches_random))
# pa_random_data <- cbind(
#   pa_random_occurrence,
#   pa_random_coords
# )

write.csv(
  x = pa_random_data,
  file = "data/tabular/presence_absence_random_sampling.csv",
  row.names = FALSE
)

# presence-absence data with sampling locations biased towards areas closer to
# major cities
pa_bias_data <- as_tibble(catches_bias_weighted) %>%
  dplyr::select(count, presence) %>%
  bind_cols(as_tibble(crds(catches_bias_weighted)))

write.csv(
  x = pa_bias_data,
  file = "data/tabular/presence_absence_bias_sampling.csv",
  row.names = FALSE
)

# presence-absence data with sampling locations biased towards areas with higher
# abundance
pa_bias_abund_data <- as_tibble(catches_abundance_weighted) %>%
  dplyr::select(count, presence) %>%
  bind_cols(as_tibble(crds(catches_abundance_weighted)))

write.csv(
  x = pa_bias_abund_data,
  file = "data/tabular/presence_absence_bias_abund_sampling.csv",
  row.names = FALSE
)

# to do:

# output 4x rasters:
#  rel_abund
#  prob_present,
#  bias
#  reported_occurrence_rate

# output 4x datasets:
#  occurrence_coords (presence-only coordinates biased towards major cities)
#  pa_random_data (presence-absence, locations randomly sampled)
#  pa_bias_data (presence-absence, locations biased towards major cities)
#  pa_bias_abund_data (presence-absence, locations biased towards higher
#      abundance areas)
