library(terra)
library(tidyverse)

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

rastpointplot <- function(
  r,
  v,
  pch = 16,
  cex = 0.5
){
  
  plot(r)
  points(v, pch = pch, cex = cex)
  
}

