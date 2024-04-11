# preparing data part - simulating joint PA and PO data for multiple species
# with shared bias and predictive covariates

# load packages
library(tidyverse)
library(terra)
source("R/functions.R")

# read our rasters in
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

writeRaster(covs,
            "data/grids/covs_mspp.tif",
            overwrite = TRUE)

# bias layer
bias <- rescale_travel ^ 2
names(bias) <- "bias"

writeRaster(bias,
            "data/grids/bias_mspp.tif",
            overwrite = TRUE)

# generate the fake 'true' relative abundance raster for the target species

# generate an unscaled relative abundance
# rel_abund_unscaled <- exp(-1 + covs$tmax * 0.1)

# create your own here:

rel_abund_unscaled <- exp(-1 +
                            covs$tseas *  -0.03   +
                            covs$tmax *  -0.02   +
                            covs$trange *  0.01   +
                            covs$trange ^ 2 *  0  )

# rescale the relative abundance, from 0 to 1
rel_abund <- rescale_abundance(rel_abund_unscaled)
names(rel_abund) <- "relative_abundance"

# sample abundance data at a random set of locations
n_samples <- 100

# random locations all over the country - unweighted sampling
sample_locations_random <- random_locations(kenya_mask,
                                            n_samples,
                                            weighted = FALSE)

# simulate PA data
catches_random <- sim_catches(sample_locations = sample_locations_random,
                              relative_abundance = rel_abund,
                              max_average_catch_size = 10)


# simulate PO data
prob_present <- probability_of_presence(rel_abund,
                                        max_average_catch_size = 10)
names(prob_present) <- "prob_present"

reported_occurrence_rate <- prob_present * bias
names(reported_occurrence_rate) <- "reported_occurrence_rate"

n_occurrences <- 300
sample_locations_bias_weighted <- random_locations(reported_occurrence_rate,
                                                   n_occurrences)

pa_random_data <- as_tibble(catches_random) %>%
  dplyr::select(count, presence) %>%
  bind_cols(as_tibble(crds(catches_random))) %>%
  mutate(
    site_id = row_number(),
    .before = everything()
  )

dir.create("data/tabular", showWarnings = FALSE)
write.csv(
  x = pa_random_data,
  file = "data/tabular/presence_absence_random_sampling.csv",
  row.names = FALSE
)

write.csv(
  x = crds(sample_locations_bias_weighted),
  file = "data/tabular/presence_only_points.csv",
  row.names = FALSE
)


writeRaster(rel_abund,
            "data/grids/rel_abund_mspp.tif",
            overwrite = TRUE)

writeRaster(prob_present,
            "data/grids/prob_present_mspp.tif",
            overwrite = TRUE)

writeRaster(reported_occurrence_rate,
            "data/grids/reported_occurrence_rate_mspp.tif",
            overwrite = TRUE)


# add fake other species data - all with restricted distributions (different
# covariates) and the same bias as the target species

n_species_each <- 10
focal_species_n_points <- 10 + rpois(n_species_each, rlnorm(n_species_each, 5, 1))

# random presence_absence locations
# pa_random_data <- read_csv("data/tabular/presence_absence_random_sampling.csv")
pa_coords <- pa_random_data[, c("x", "y")]

# now do focal ones, each time selecting a different bioclim layer and including it in the sampling

focal_occ_df <- tibble(x = numeric(0),
                       y = numeric(0),
                       species_id = integer(0))

focal_pa_df <- tibble(site_id = integer(0),
                      presence = integer(0),
                      species_id = integer(0))

for(i in seq_len(n_species_each)) {
  
  which_cov <- sample.int(dim(covs)[3], 1)
  prob_pres <- app(rnorm(1) * scale(covs[[which_cov]]), plogis)
  rep_occ_rate <- prob_pres * bias
  
  # simulate presence ony data
  one_focal_species_occ_points <- terra::spatSample(rep_occ_rate,
                                                    focal_species_n_points[i],
                                                    method = "weights",
                                                    replace = TRUE,
                                                    na.rm = TRUE,
                                                    as.points = TRUE) 
  
  one_focal_species_df <- one_focal_species_occ_points %>%
    crds() %>%
    as_tibble() %>%
    mutate(
      species_id = letters[i]
    )
  
  focal_occ_df <- rbind(focal_occ_df, one_focal_species_df)
  
  # simulate presence-absence data
  p <- terra::extract(prob_pres, pa_coords)[, 2]
  presence <- rbinom(length(p), 1, p)
  
  focal_pa_df <- rbind(focal_pa_df,
                       tibble(
                         site_id = seq_along(p),
                         presence = presence,
                         species_id = letters[i])
  )
  
}

focal_pa_tabular <- focal_pa_df %>%
  pivot_wider(names_from = species_id,
              values_from = presence) %>%
  left_join(
    bind_cols(site_id = seq_len(nrow(pa_coords)),
              pa_coords),
    
    by = "site_id"
  ) %>%
  relocate(x, y,
           .after = site_id)

# save these items
write.csv(focal_occ_df,
          file = "data/tabular/other_species_po_data.csv",
          row.names = FALSE)

write.csv(focal_pa_tabular,
          file = "data/tabular/other_species_pa_data.csv",
          row.names = FALSE)
