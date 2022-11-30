# first models
library(terra)
library(tidyverse)
library(predicts) # remotes::install_github("rspatial/predicts")

# sample random background points
n_background_points <- 100

random_bg <- terra::spatSample(
  x = kenya_mask,
  size = n_background_points,
  na.rm = TRUE,
  as.points = TRUE
)

# our presence data is collected in biased fashion
sample_locations_bias_weighted

#let's just use the 3 key covariates for now:
covs
plot(covs)


# combine data for analysis
sdm_data_1 <- make_sdm_data(
  presences = sample_locations_bias_weighted,
  absences = random_bg,
  covariates = covs
)


# fit a model!

model_1 <- glm(
  occ ~ tseas + tmax + trange,
  data = sdm_data_1,
  family = binomial()
)

summary(model_1)

# look at the partial_response

pr1 <- partialRe
