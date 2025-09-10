# count model
model_count_random_logistic <- glm(
  count ~ tseas + tmax + trange,
  data = data_pa_random,
  family = stats::poisson()
)
summary(model_count_random_logistic)

# writing out the glm above
library(greta)
# data
tseas <- data_pa_random$tseas
tmax <- data_pa_random$tmax
trange <- data_pa_random$trange
presence <- data_pa_random$presence

# covariates
beta_tseas <- normal(0, 1)
beta_tmax <- normal(0, 1)
beta_trange <- normal(0, 1)
# intercept
alpha <- normal(0,1)

# link scale probability
gamma <- alpha + beta_tseas*tseas + beta_tmax*tmax + beta_trange*trange

logitfun <- function(gamma){
  1/(1+exp(-gamma))
}

# response scale probability
prob <- logitfun(gamma)

# likelihood
presence ~ binomial(1, prob = prob)

# model
m <- model(beta_tseas, beta_tmax, beta_trange)

draws <- mcmc(m, n_samples = 500)

library(bayesplot)
summary(draws)
mcmc_intervals(draws)


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