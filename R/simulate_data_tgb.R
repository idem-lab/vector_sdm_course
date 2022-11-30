# preparing data part 2

# should run without part 1 if data/grids/ contains the below rasters.
# as an alternative to runnning script preparing_data_1.R, grids can be
# downloaded into data/grids from:
# https://melbourne.figshare.com/articles/dataset/vector_sdm_course_grids/
# you can do this manually or work out how to in R!

# load packages
library(tidyverse)
library(terra)

# read our rasters in
kenya_mask <- terra::rast("data/grids/kenya_mask.tif")
bc_kenya <- terra::rast("data/grids/bc_kenya.tif")
rescale_travel <- terra::rast("data/grids/rescale_travel.tif")


# add fake target species data - make some that are widespread, some that are
# more restricted (add that information as a column). We can be confident that
# all have the same sampling bias as the target species (fake species made by
# the user)

# sample only biased ones from rescale_travel

# sample biased and environmentally specific ones from rescale_travel
# and some bioclim variables

# make two sets of fake species, with different numbers of points
n_species_each <- 10
widespread_species_n_points <- 10 + rpois(n_species_each, rlnorm(n_species_each, 5, 1))
focal_species_n_points <- 10 + rpois(n_species_each, rlnorm(n_species_each, 5, 1))

# simulate all the widespread ones just from the bias layer
widespread_all_points <- terra::spatSample(rescale_travel,
                                           sum(widespread_species_n_points),
                                           method = "weights",
                                           replace = TRUE,
                                           na.rm = TRUE,
                                           as.points = TRUE)

# plot(rescale_travel)
# points(widespread_all_points, cex = 0.2)

# label them randomly acorrding to different species
widespread_df <- widespread_all_points %>%
  crds() %>%
  as_tibble() %>%
  mutate(
    species_id = rep(seq_len(n_species_each),
                     widespread_species_n_points),
    type = "widespread"
  )


# now do focal ones, each time selecting a different bioclim layer and including it in the sampling

focal_df <- tibble(x = numeric(0),
                   y = numeric(0),
                   species_id = integer(0),
                   type = character(0))

for(i in seq_len(n_species_each)) {
  
  which_bioclim <- sample.int(dim(bioclim_kenya)[3], 1)
  layer <- bioclim_kenya[[which_bioclim]] * rescale_travel
  
  layer[] <- layer[] - min(layer[], na.rm = TRUE)
  layer[] <- layer[]/ max(layer[], na.rm = TRUE)
  
  # simulate all the widespread ones just from the bias layer
  one_focal_species_points <- terra::spatSample(layer ^ 2,
                                                focal_species_n_points[i],
                                                method = "weights",
                                                replace = TRUE,
                                                na.rm = TRUE,
                                                as.points = TRUE) 
  
  one_focal_species_df <- one_focal_species_points %>%
    crds() %>%
    as_tibble() %>%
    mutate(
      species_id = i,
      type = "focal"
    )
  
  focal_df <- rbind(focal_df, one_focal_species_df)
  
}

species_df <- rbind(focal_df, widespread_df)



# save these items
save(
  list = c(
    "widespread_df",
    "focal_df",
    "species_df"
  ),
  file = "data/point_dfs.RData"
)
