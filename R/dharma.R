# dharma
library(DHARMa)


simulationOutput <- simulateResiduals(fittedModel = m, plot = FALSE)

# plot simulated residuals
plot(simulationOutput)
# err looks really friggin good

# check if resids look normal
hist(residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7)))
# yep

# test overdispersion
testDispersion(simulationOutput)
# look way overdispersed, test non sig - low power?

# test for zero inflation
testZeroInflation(simulationOutput)
# nah looks fine as - thanks NB?
