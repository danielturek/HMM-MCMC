## Latent State MCMC for the dipper model

## load nimble library
library(nimble)

## load data
load('../data/dipperData.RData')

## define hierarchical model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for (i in 1:nind) {
        x[i, first[i]] <- 1
        for (t in (first[i] + 1):k) {
            mu_x[i, t] <- phi * x[i, t-1]
            mu_y[i, t] <- p * x[i, t]
            x[i, t] ~ dbin(mu_x[i, t], 1)
            y[i, t] ~ dbin(mu_y[i, t], 1)
        }
    }
})

## define model constants, data, and initial values
constants <- list(k=k, nind=nind, first=first)
data <- list(y=y)
inits <- list(phi=0.6, p=0.9, x=x_init)

## create model object
Rmodel <- nimbleModel(code, constants, data, inits)

## specify MCMC algorithm
spec <- configureMCMC(Rmodel)
spec$printSamplers()

## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC
set.seed(0)
Cmcmc$run(10000)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

