library(tidyverse)
library(terra)
#library(malariaAtlas)
#library(sf)
library(geodata)



# get extent map for our basic analyses - using outline of Kenya

kenya_vector <- gadm(
  country = "KEN",
  level = 0,
  path = "data/"
) # slow to run though not large

kenya_vector

plot(kenya_vector)

# let's download some data

# travel time from city to cell
# you can explore the sizes of city using ?travel time, and change size if you wish
travel <- travel_time(
  to = "city",
  size = 2,
  path = "data/"
) # 421 MB - may be slow to download

# it's a good idea to have a look at your objects and even plot them


# bioclimactic variables from worldclim
# https://worldclim.org/data/bioclim.html
bioclim_kenya <- worldclim_country(
  country = "KEN",
  var = "bio",
  res = 0.5,
  path = "data/"
)

# have a look at this object first
# Do you think you don't want to plot all layers at once?

# process these data to match our extent

# generate mask of our area of interest so we can use it to process other data into this shape
kenya_mask <- bioclim_kenya[[1]] %>%
  mask(kenya_vector) * 0 + 1

plot(kenya_mask)


# process bioclim into shape
bc_kenya <- mask(
  x = bioclim_kenya,
  mask = kenya_mask
)

# process travel data
# this is a very large raster, we speed up the operation by cropping then masking
# try looking at the result of only cropping to understand why both steps are needed
travel_kenya <- crop(
  x = travel,
  y = kenya_mask
) %>%
  mask(
    mask = kenya_mask
  )

# let's have a look at it


# now let's create a bias layer out of our travel data
# it should scale from 0-1 where 0 is hard to get to and 1 is easy

plot(travel_kenya)

rescale_travel <- travel_kenya

rescale_travel[] <- rescale_travel[]/ max(rescale_travel[], na.rm = TRUE)

rescale_travel[] <- 1 - rescale_travel[]

plot(rescale_travel)

# save these outputs
save(
  list = c(
    "kenya_mask",
    "bc_kenya",
    "rescale_travel"
  ),
  file = "data/basic_data.RData"
) # check that this actually saves properly (i.e. does it have the in-memory issue that raster:: does)


#plot(bc_kenya[[1:10]])
#plot(bc_kenya[[11:19]])
# 1, 5, 6, 8, 9, 10, 11 + 1; 12, 13, 14, 16, 17, 18, 19 -1
# 2, 7, + 1
# 3, 4, -1

# BIO1 = Annual Mean Temperature
# 
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# 
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# 
# BIO4 = Temperature Seasonality (standard deviation ×100)
# 
# BIO5 = Max Temperature of Warmest Month
# 
# BIO6 = Min Temperature of Coldest Month
# 
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# 
# BIO8 = Mean Temperature of Wettest Quarter
# 
# BIO9 = Mean Temperature of Driest Quarter
# 
# BIO10 = Mean Temperature of Warmest Quarter
# 
# BIO11 = Mean Temperature of Coldest Quarter
# 
# BIO12 = Annual Precipitation
# 
# BIO13 = Precipitation of Wettest Month
# 
# BIO14 = Precipitation of Driest Month
# 
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# 
# BIO16 = Precipitation of Wettest Quarter
# 
# BIO17 = Precipitation of Driest Quarter
# 
# BIO18 = Precipitation of Warmest Quarter
# 
# BIO19 = Precipitation of Coldest Quarter






# add fake target species data - make some that are widespread, some that are
# more restricted (add that information as a column). We can be confident that
# all have the same sampling bias as the target species (fake species made by
# the user)

# sample only biased ones from rescale_travel

# sample biased and environmentally specific ones from rescale_travel and some bioclim variables

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



  
  
  