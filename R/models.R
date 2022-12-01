# fit models

library(terra)
library(tidyverse)
#library(predicts) # remotes::install_github("rspatial/predicts")
library(glmnet)
library(maxnet)
source("R/functions.R")

# read in our data!

# rasters
rel_abund <- terra::rast("data/grids/rel_abund.tif")
prob_present <- terra::rast("data/grids/prob_present.tif")
bias <- terra::rast("data/grids/bias.tif")
reported_occurrence_rate <- terra::rast("data/grids/reported_occurrence_rate.tif")

kenya_mask <- terra::rast("data/grids/kenya_mask.tif")
bc_kenya <- terra::rast("data/grids/bc_kenya.tif")
rescale_travel <- terra::rast("data/grids/rescale_travel.tif")

# points 
occurrence_coords <- read_csv("data/tabular/presence_only_data.csv")
pa_random_data <- read_csv("data/tabular/presence_absence_random_sampling.csv")
pa_bias_data <- read_csv("data/tabular/presence_absence_bias_sampling.csv")
pa_bias_abund_data  <- read_csv("data/tabular/presence_absence_bias_abund_sampling.csv")

species_df <- read_csv("data/tabular/species_df.csv")

# subset bc_keyna to our covariate set (and save this again for ease of use)
covs <- bc_kenya[[c(4,5,7)]]
names(covs) <- c("tseas", "tmax", "trange")

terra::writeRaster(
  x = covs,
  filename = "data/grids/covariates.tif"
)


### Model: logistic regression of random presence-absence data

pa_random_data

data_pa_random <- model_data_presence_absence(
  pa_data = pa_random_data,
  covariates = covs
)

data_pa_random

# fit a simple model!
model_pa_random_logistic <- glm(
  presence ~ tseas + tmax + trange,
  data = data_pa_random,
  family = binomial()
)
summary(model_pa_random_logistic)


# plot the partial responses for each
# predictor variable (covariate)
partial_response_plot(
  model = model_pa_random_logistic,
  data = data_pa_random,
  var = "tmax"
)
# now do
# tseas
# trange
partial_response_plot(
  model = model_pa_random_logistic,
  data = data_pa_random,
  var = "tseas"
)

partial_response_plot(
  model = model_pa_random_logistic,
  data = data_pa_random,
  var = "trange"
)

# predict our distribution based on our model
# and covariates
pred_pa_random_logistic <- sdm_predict(
  model = model_pa_random_logistic,
  covariates = covs
)

# plot it
plot(pred_pa_random_logistic)

#compare it with the truth
plot(c(rel_abund, prob_present, pred_pa_random_logistic))
plot(c(prob_present, pred_pa_random_logistic))
### Model: logistic regression of presence-only data
# with random background 

# sample random background points
n_background_points <- 100

random_bg <- terra::spatSample(
  x = kenya_mask,
  size = n_background_points,
  na.rm = TRUE,
  as.points = TRUE
)

rastpointplot(kenya_mask, random_bg)

# put presence and background data together
# with covariates

data_po_random_bg <- model_data_presence_only(
  presences = occurrence_coords,
  absences = random_bg,
  covariates = covs
)


# fit a simple model!
model_po_random_bg_logistic <- glm(
  presence ~ tseas + tmax + trange,
  data = data_po_random_bg,
  family = binomial()
)
summary(model_po_random_bg_logistic)

# partial response of each variable
partial_response_plot(
  model = model_po_random_bg_logistic,
  data = data_po_random_bg,
  var = "tmax"
)
# do others!
partial_response_plot(
  model = model_po_random_bg_logistic,
  data = data_po_random_bg,
  var = "tseas"
)

partial_response_plot(
  model = model_po_random_bg_logistic,
  data = data_po_random_bg,
  var = "trange"
)
# predict our distribution based on our model and covariates
pred_po_random_bg_logistic <- sdm_predict(
  model = model_po_random_bg_logistic,
  covariates = covs
)

# plot it
plot(pred_po_random_bg_logistic)

plot(c(pred_pa_random_logistic, pred_po_random_bg_logistic))
# now compare that prediction with the truth
plot(c(prob_present, pred_pa_random_logistic, pred_po_random_bg_logistic, reported_occurrence_rate))


### Model:presence-only with maxnet (R version of maxent)

#use the presence only and random background again
data_po_random_bg

# model in maxnet

model_po_random_bg_maxent <- maxnet(
  p = data_po_random_bg$presence,
  data = data_po_random_bg %>% select(-presence),
  addsamplestobackground = FALSE # already created them
)


# partial response of each variable
# different for maxnet than the glms
plot(model_po_random_bg_maxent, "tseas")

plot(model_po_random_bg_maxent, "tmax")

plot(model_po_random_bg_maxent, "trange")
# do others!


# predict our distribution based on our model and covariates
pred_po_random_bg_maxent <- sdm_predict(
  model = model_po_random_bg_maxent,
  covariates = covs
)

# plot it
plot(pred_po_random_bg_maxent)

plot(c(pred_po_random_bg_logistic, pred_po_random_bg_maxent))

plot(c(prob_present, pred_pa_random_logistic,
       pred_po_random_bg_logistic, pred_po_random_bg_maxent))


### Model:presence-only with maxnet (R version of maxent) with bias layer

#use the presence only and random background again
data_po_random_bg

# extract bias values from presence and backgroud locations
maxent_bias_df <-  model_data_presence_only(
  presences = occurrence_coords,
  absences = random_bg,
  covariates = bias
)
# log them and create vector
maxent_bias <- log(maxent_bias_df$bias)

# model in maxnet

model_po_random_bg_maxent_bias <- maxnet(
  p = data_po_random_bg$presence,
  data = data_po_random_bg %>% select(-presence),
  offset = maxent_bias %>% as.matrix(),
  addsamplestobackground = FALSE # becase we have included background
)


# partial response of each variable
# different for maxnet than the glms
plot(model_po_random_bg_maxent_bias, "tseas")
# do others!


# predict our distribution based on our model and covariates
pred_po_random_bg_maxent_bias <- sdm_predict(
  model = model_po_random_bg_maxent_bias,
  covariates = covs
)

# plot it
plot(pred_po_random_bg_maxent_bias)


### Model: glm with target group-background

# our data
# presences
occurrence_coords

# this time we will use other species as "absences"
species_df



data_po_tgb_all <- model_data_presence_only(
  presences = occurrence_coords,
  absences = species_df %>% dplyr::select(x, y),
  covariates = covs
)


# fit a simple model!
model_po_tgb_logistic <- glm(
  presence ~ tseas + tmax + trange,
  data = data_po_tgb_all,
  family = binomial()
)
summary(model_po_tgb_logistic)

# partial response of each variable
partial_response_plot(
  model = model_po_tgb_logistic,
  data = data_po_tgb_all,
  var = "tmax"
)

# predict our distribution based on our model and covariates
pred_po_tgb_logistic <- sdm_predict(
  model = model_po_tgb_logistic,
  covariates = covs
)

# plot it
plot(pred_po_tgb_logistic)

# now compare that prediction with the truth
plot(c(pred_pa_random_logistic, pred_po_tgb_logistic))
