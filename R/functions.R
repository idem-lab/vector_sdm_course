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
  min_value <- global(unscaled_abundance, "min", na.rm = TRUE)[1, 1]
  unscaled_abundance <- unscaled_abundance - min_value + 0.00001
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

# plot a raster with simple points over it
rastpointplot <- function(
  r,
  v,
  pch = 16,
  cex = 0.5
){
  
  plot(r)
  points(v, pch = pch, cex = cex)
  
}


# # combine sdm data from spatial to df
# make_sdm_data <- function(
#   presences,
#   absences,
#   covariates
# ){
#   
#   pvals <- terra::extract(covariates, presences)
#   avals <- terra::extract(covariates, absences)
#   
#   rbind(
#     pvals %>%
#       as_tibble %>%
#       dplyr::select(-ID) %>%
#       mutate(occ = 1),
#     avals %>%
#       as_tibble %>%
#       dplyr::select(-ID) %>%
#       mutate(occ = 0)
#   ) %>%
#     as_tibble
#   
# }


# name prediction layer 
sdm_predict <- function(
  model,
  covariates,
  type = NULL,
  layer_name = "predicted_distribution"
){
  
  if(is.null(type)){
    if (inherits(model, "maxnet")) {
      type <- "logistic"
    } else {
      type <- "response"
    }
  }
  
  prediction <- predict(covariates, model, na.rm = TRUE, type = type)
  names(prediction) <- layer_name
  
  return(prediction)
}

# function to take presence only and background data
# and covariate rasters and form into a data frame 
# for modelling
model_data_presence_only <- function(
  presences,
  absences,
  covariates
){
  
  pvals <- terra::extract(covariates, presences)
  avals <- terra::extract(covariates, absences)
  
  rbind(
    pvals %>%
      as_tibble %>%
      dplyr::select(-ID) %>%
      mutate(presence = 1),
    avals %>%
      as_tibble %>%
      dplyr::select(-ID) %>%
      mutate(presence = 0)
  ) %>%
    as_tibble
  
}

# function to take presence absence data
# and covariate rasters and form into a data frame 
# for modelling

model_data_presence_absence <- function(
  pa_data,
  covariates
){
  
  vals <- terra::extract(
    covariates,
    pa_data %>%
      dplyr::select(x, y)
  )
  
  cbind(
    pa_data %>%
      dplyr::select(-x, -y),
    vals %>%
      as_tibble %>%
      dplyr::select(-ID)
  ) %>%
    as_tibble
  
}



# fit a partial response curve
partial_response <- function (model, data, var, type = c("response", "link"), rng = NULL, nsteps = 25) {
  
  type <- match.arg(type)
  if (missing(var)) {
    var <- names(data)[1]
  }
  else if (is.numeric(var)) {
    stopifnot(var > 0 & var <= ncol(data))
    var <- names(data)[var]
  }
  else {
    stopifnot(var %in% names(data))
  }
  if (is.factor(data[[var]])) {
    steps <- levels(data[[var]])
  }
  else {
    if (is.null(rng)) {
      rng <- range(data[[var]])
    }
    increment <- (rng[2] - rng[1])/(nsteps - 2)
    steps <- seq(rng[1] - increment, rng[2] + increment, 
                 increment)
  }
  res <- rep(NA, length(steps))
  for (i in 1:length(steps)) {
    data[[var]] <- steps[i]
    p <- predict(model, data, type = type)
    res[i] <- mean(p)
  }
  x <- data.frame(steps, res)
  names(x) <- c("var", "p")
  x
}

# plot partial response curve
partial_response_plot <- function(
  model,
  data,
  var
){
  plot(
    partial_response(
      model = model,
      data = data,
      var = var
    ),
    type = "l"
  )
}



