

predictor_values <- as.matrix(covs) %>%
  as_tibble %>%
  filter(!is.nan(tmax))

predictor_values <- predictor_values[floor(seq(from = 1, to = nrow(predictor_values), length.out = 20)),]

predictor_values


coefs <- model_pa_random_logistic$coefficients %>%
  as.list

coefs$intercept <- coefs$`(Intercept)`

predictor_values %>%
  mutate(
    # model formila presence ~ tseas + tmax + trange
    prediction_link = coefs$intercept +
      tseas * coefs$tseas +
      tmax * coefs$tmax +
      trange * coefs$trange,
    prediction_response = 1/(1+exp(-prediction_link))
  )



# calculate by hand

# our data
tseas <- data_pa_random$tseas
tmax <- data_pa_random$tmax
trange <- data_pa_random$trange
#presence <- data_pa_random$presence

# some coefficients/ slopes - play around with these
beta_tseas <- 0.001
beta_tmax <- 0.002
beta_trange <- 0.1
alpha <- -4

# calculate probability on link scale
gamma <- alpha + beta_tseas*tseas + beta_tmax*tmax + beta_trange*trange

# logit link function
logitfun <- function(gamma){
  1/(1+exp(-gamma))
}

# calculate probability of presence
logitfun(gamma)

# calculate a partial response curve
tseas_mean <- mean(data_pa_random$tseas)
tmax_mean <- mean(data_pa_random$tmax)
trange_mean <- mean(data_pa_random$trange)

pr_tmax_data <- seq(
  from = 20,
  to = 40,
  by = 0.2
)

# some coefficients/ slopes - play around with these
beta_tseas <- 0.001
beta_tmax <- 0.2
beta_trange <- 0.1
alpha <- -5

# calculate probability on link scale
gamma_pr <- alpha + beta_tseas*tseas_mean + beta_tmax*pr_tmax_data + beta_trange*trange_mean

# logit link function
logitfun <- function(gamma){
  1/(1+exp(-gamma))
}

# calculate probability of presence
prob_pr <- logitfun(gamma_pr)

plot(
  x = pr_tmax_data,
  y = prob_pr,
  ylim = c(0, 1)
)
