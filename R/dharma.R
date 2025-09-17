# dharma
library(DHARMa)

## well fitting model

# assess the residual
simulationOutput_po_rbg_logistic <- simulateResiduals(
  fittedModel = model_pa_random_logistic,
  plot = FALSE
)

# plot simulated residuals
plot(simulationOutput_po_rbg_logistic)

# check if resids look normal
simulationOutput_po_rbg_logistic |>
  residuals(
    quantileFunction = qnorm,
    outlierValues = c(-7,7)
  ) |>
  hist()


# test overdispersion
testDispersion(simulationOutput_po_rbg_logistic)

# test for zero inflation
testZeroInflation(simulationOutput_po_rbg_logistic)


## Poorly specified model
summary(model_pa_random_poisson_quad)

# assess the residual

# assigning the model per models.R causes R to abort
# simulationOutput_pa_random_poisson_quad <- simulateResiduals(
#   fittedModel = model_pa_random_poisson_quad,
#   plot = FALSE
# )

# this is a hack to get around this issue whil still producing the
# awful model
simulationOutput_pa_random_poisson_quad <- simulateResiduals(
  fittedModel = glm(
    presence ~ I(tmax^2),
    data = data_pa_random,
    family = stats::poisson()
  ),
  plot = FALSE
)

# plot simulated residuals
plot(simulationOutput_pa_random_poisson_quad)

# check if resids look normal
simulationOutput_pa_random_poisson_quad |>
  residuals(
    quantileFunction = qnorm,
    outlierValues = c(-7,7)
  ) |>
  hist()


# test overdispersion
testDispersion(simulationOutput_pa_random_poisson_quad)

# test for zero inflation
testZeroInflation(simulationOutput_pa_random_poisson_quad)

