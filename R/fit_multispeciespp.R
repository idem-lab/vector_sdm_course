# fit a multispecies point process model from Fithian et al.

# install the package from github, using the 'remotes' package
# install.packages("remotes")
# library(remotes)
# install_github("wfithian/multispeciesPP")

library(terra)
library(tidyverse)
library(multispeciesPP)

source("R/functions.R")

# load rasters for fitting
kenya_mask <- terra::rast("data/grids/kenya_mask.tif")
rescale_travel <- terra::rast("data/grids/rescale_travel.tif")
covs <- terra::rast("data/grids/covs_mspp.tif")
ext(covs) <- ext(kenya_mask)
crs(covs) <- crs(kenya_mask)


# rasters for comparison
prob_present <- terra::rast("data/grids/prob_present_mspp.tif")
names(prob_present) <- "prob_present"
ext(prob_present) <- ext(kenya_mask)
rep_occ_rate <- terra::rast("data/grids/reported_occurrence_rate_mspp.tif")
names(rep_occ_rate) <- "rep_occ_rate"
ext(rep_occ_rate) <- ext(kenya_mask)

# load the presence-absence and presence-only data for the target and other species

# for the target species
pa_target_species <- read.csv("data/tabular/presence_absence_random_sampling.csv")
po_target_species <- read.csv("data/tabular/presence_only_points.csv")

# for other species with the same bias
pa_other_species <- read.csv("data/tabular/other_species_pa_data.csv")
po_other_species <- read.csv("data/tabular/other_species_po_data.csv")

# generate random locations for background data
n_background_points <- 1000
random_bg_coords <- terra::spatSample(
  x = kenya_mask,
  size = n_background_points,
  na.rm = TRUE,
  as.points = TRUE
) %>%
  crds()


# fit a presence-only random BG model for just the target species
data_po_random_bg <- model_data_presence_only(
  presences = as_tibble(po_target_species),
  absences = as_tibble(random_bg_coords),
  covariates = covs
)

# fit a simple model!
model_po_random_bg_logistic <- glm(
  presence ~ tseas + tmax + trange,
  data = data_po_random_bg,
  family = binomial()
)

# predict our distribution based on our model and covariates
pred_po_random_bg_logistic <- sdm_predict(
  model = model_po_random_bg_logistic,
  covariates = covs
)

plot(pred_po_random_bg_logistic)


# format these data into the format expected by multispeciesPP
pa_coords <- pa_target_species %>%
  select(x, y)
pa_covariates <- terra::extract(covs, pa_coords) %>%
  select(-ID)

# presence-absence data should be a single dataframe with presence/absence for
# all the species at the same set of sites, and all the covariates added
pa_data <- pa_target_species %>%
  # we don't use the count data
  select(-count, -x, -y) %>%
  # the presence colummn should be our species' name
  rename(target = presence) %>%
  # combine with the other species
  left_join(
    pa_other_species,
    by = "site_id"
  ) %>%
  # we don't need site_id or covariates any more
  select(-site_id, -x, -y) %>%
  relocate(
    target,
    .before = "a"
  ) %>%
  # add on the covariate values
  bind_cols(
    pa_covariates, .)

# presence-only data should be a named list, with the species names, containing
# the covariate values for both the distribution and bias models at the
# occurrence locations

#  first extract the covariate and bias values for target and other species

covs_and_bias <- c(covs, rescale_travel)

po_all_species <- bind_rows(
  po_target_species %>%
    mutate(species_id = "target"),
  po_other_species
)

po_covs_all_species <- po_all_species %>%
  select(x, y) %>%
  terra::extract(covs_and_bias, .) %>%
  select(-ID) %>%
  mutate(species_id = po_all_species$species_id)

# convert to a named list
po_covs_all_species_list <- po_covs_all_species %>%
  group_by(species_id) %>%
  summarise(named_vec = list(.)) %>%
  deframe() %>%
  lapply(ungroup)

# subset each one
po_covs_all_species_list <- mapply(
  function(name, tibble) {
    filter(tibble, species_id == name)
  },
  names(po_covs_all_species_list),
  po_covs_all_species_list,
  SIMPLIFY = FALSE
)

# background coordinate values
bg <- terra::extract(covs_and_bias, random_bg_coords)

n_pa_obs <- 50
pa_data_keep_idx <- sample.int(nrow(pa_data), n_pa_obs)

full.mod <- multispeciesPP(sdm.formula = ~ tseas + tmax + trange,
                           bias.formula = ~travel_time_to_cities_2,
                           PA = pa_data[pa_data_keep_idx, ],
                           PO = po_data_list,
                           BG = bg)

# prediction is more difficult
non_na_idx <- which(!is.na(as.vector(kenya_mask)))
pred_vals <- covs_and_bias[non_na_idx]

link_preds <- multispeciesPP::predict.multispeciesPP(full.mod,
                                                     newdata = pred_vals,
                                                     species = "target")
prob_preds <- 1 - exp(-exp(link_preds))

predicted_distribution <- kenya_mask
names(predicted_distribution) <- "predicted_distribution"
predicted_distribution[non_na_idx] <- prob_preds[, 1]

# plot(predicted_distribution)

plot(c(prob_present,
       rep_occ_rate,
       predicted_distribution,
       pred_po_random_bg_logistic))

