library(terra)

# Given a raster and a number of sampling locations to simulate, return a
# SpatVector object of the locations. If weighted = TRUE, then treat the values
# of the raster as the relative number of points to put in each cell. If replace
# = TRUE, then multiple samples can appear in the same raster cell.
random_locations <- function(raster, n_samples, weighted = TRUE, replace = TRUE) {
  
  if (weighted) {
    method <- "weights"
  } else {
    method <- "random"
  }
  
  terra::spatSample(raster,
                    n_samples,
                    method = method,
                    replace = replace,
                    na.rm = TRUE,
                    as.points = TRUE)
  
}

# given a SpatVector object of locations at which to sample, and a SpatRaster
# object of the relative abundance (must be scaled to a maximum value of 1), and
# a maximum average catch size (the average catch size at the location with the
# highest abundance), return the SpatVector object with extra columns for the a
# simulated catch size 'count', and binary variable for whether the species was
# present in the catch
sim_catches <- function(sample_locations,
                        relative_abundance,
                        max_average_catch_size = 5000) {
  # how many samples?
  n_samples <- nrow(sample_locations)
  
  # if you were to trap here repeatedly, what would the average catch size be,
  # in the long term?
  average_catch_size <- relative_abundance * max_average_catch_size
  
  # get the values for this at the sample locations
  expected_catch_size <- terra::extract(average_catch_size,
                                        sample_locations)[, 2]
  
  # given this average catch size, sample a single catch (assuming a Poisson
  # distribution for catch sizes)
  sample_locations$count <- rpois(n_samples, expected_catch_size)
  
  # record a 1 if any were caught, and a 0 otherwise
  sample_locations$presence <- pmin(sample_locations$count, 1)
  
  # return the locations with this info attached
  sample_locations
  
}

# given an unscaled relative abundance raster, scale it to have maximum value of 1 
rescale_abundance <- function(unscaled_abundance) {
  max_value <- global(unscaled_abundance, "max", na.rm = TRUE)[1, 1]
  unscaled_abundance / max_value
}

# given a SpatRaster object of the relative abundance (must be scaled to a
# maximum value of 1), and a maximum average catch size (the average catch size
# at the location with the highest abundance), return a SpatRaster of the
# probability of the species being present in a random catch at each location
probability_of_presence <- function(relative_abundance,
                                    max_average_catch_size = 5000) {
  
  # if you were to trap here repeatedly, what would the average catch size be,
  # in the long term?
  average_catch_size <- relative_abundance * max_average_catch_size
  
  # what is the probability of detecting one or more mosqiutoes in a poisson
  # sample with an average catch of this size
  probability_one_or_more <- 1 - exp(-(average_catch_size))
  
  probability_one_or_more
  
}

# fake rasters to get the code working - these to be replaced with real ones

# fake raster from the terra package
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

# mask layer (useful for processing)
mask <- r * 0

# make a fake bias layer
bias <- mask
cell_ids <- cells(mask)
bias[cell_ids] <- terra::xFromCell(mask, cell_ids)
bias <- mask(bias, mask)
bias <- bias - global(bias, "min", na.rm = TRUE)[1,1]
bias <- bias / global(bias, "max", na.rm = TRUE)[1,1]
bias <- bias ^ 2

# make some fake covariates
cov_a <- r
cov_b <- -log(r)
cov_c <- app(r / global(r, "max", na.rm = TRUE)[1, 1], qlogis)
covs <- c(cov_a, cov_b, cov_c)
names(covs) <- c("a", "b", "c")


# generate the fake 'true' relative abundance raster

# generate an unscaled relative abundance

rel_abund_unscaled <- exp(-1 + covs$a * 0.1)

# create your own here:

# rel_abund_unscaled <- exp(? +
#                            covs$a *  ?   +
#                            covs$b *  ?   +
#                            covs$c *  ?   + 
#                            covs$c ^ 2 *   ?  )

# rescale the relative abundance, from 0 to 1
rel_abund <- rescale_abundance(rel_abund_unscaled)


# sample abundance data at a random set of locations
n_samples <- 100

# random locations all over the country - unweighted sampling
sample_locations_random <- random_locations(mask,
                                            n_samples,
                                            weighted = FALSE)

catches_random <- sim_catches(sample_locations = sample_locations_random,
                              relative_abundance = rel_abund)

plot(rel_abund)
points(catches_random, pch = 21, bg = catches_random$presence)


# random locations, biased as per the bias layer - e.g. convenience samples
sample_locations_bias_weighted <- random_locations(bias,
                                                   n_samples)

catches_bias_weighted <- sim_catches(sample_locations = sample_locations_bias_weighted,
                                     relative_abundance = rel_abund)

plot(rel_abund)
points(catches_bias_weighted, pch = 21, bg = catches_bias_weighted$presence)



# random locations, biased as per the relative abundance layer - e.g. targetted
# to areas of high abundance (where malaria interventions happen?)
sample_locations_abundance_weighted <- random_locations(rel_abund ^ (1/3),
                                                        n_samples)

catches_abundance_weighted <- sim_catches(sample_locations = sample_locations_abundance_weighted,
                                          relative_abundance = rel_abund)

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
                                     relative_abundance = rel_abund)
occurrence_coords <- crds(catches_bias_weighted[catches_bias_weighted$presence == 1])

plot(mask)
points(occurrence_coords, pch = 16)

# this is great, but when simulating, it's difficult to get a realistic number
# of occurrence records for model fitting! You can try changing the number of
# sampling locations (n_samples in the code above) until it looks good:

sample_locations_bias_weighted <- random_locations(bias,
                                                   2000)  # <- change this number
catches_bias_weighted <- sim_catches(sample_locations = sample_locations_bias_weighted,
                                     relative_abundance = rel_abund)
occurrence_coords <- crds(catches_bias_weighted[catches_bias_weighted$presence == 1])

# number of records
nrow(occurrence_coords)
plot(mask)
points(occurrence_coords, pch = 16)


# 2. alternatively, you can simulate the distribution of occurrence points by
# simulating locations, biased by the *product* of the bias and probability of
# detection in each cell. The probability of detection is calculated using a
# formula that assumes the occurrence comes from a catch (like the simulations
# above), and that the number of mosquitoes caught is a Poisson sample, given
# the average number. It is the probability of observing one *or more*
# mosquitoes in the catch. 

prob_present <- probability_of_presence(rel_abund)
plot(prob_present)

reported_occurrence_rate <- bias * prob_present
n_occurrences <- 100
sample_locations_bias_weighted <- random_locations(reported_occurrence_rate,
                                                   n_occurrences)

# plot the reported occurrence rates
plot(reported_occurrence_rate)
points(sample_locations_bias_weighted, pch = 16)

