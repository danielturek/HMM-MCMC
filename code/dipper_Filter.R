## Filter MCMC for the dipper model

## load nimble library
library(nimble)

## define custom distribution
dCJS2 <- nimbleFunction(
    run = function(x = double(), length = double(), nSightingsMinus1 = double(), nNonSightings = double(), log_phi = double(), log_pp = double(), log_1minusP = double(), logChi = double(), log.p = double()) {
        logL <- length*log_phi + nSightingsMinus1*log_pp + nNonSightings*log_1minusP + logChi
        returnType(double())
        return(logL)
    }
)

rCJS2 <- nimbleFunction(
    run = function(n = integer(), length = double(), nSightingsMinus1 = double(), nNonSightings = double(), log_phi = double(), log_pp = double(), log_1minusP = double(), logChi = double()) {
        returnType(double())
        return(1)
    }
)

registerDistributions(list(
    dCJS2 = list(
        BUGSdist = 'dCJS2(length, nSightingsMinus1, nNonSightings, log_phi, log_pp, log_1minusP, logChi)',
        types = c('value = double()', 'length = double()', 'nSightingsMinus1 = double()', 'nNonSightings = double()', 'log_phi = double()', 'log_pp = double()', 'log_1minusP = double()', 'logChi = double()'),
        discrete = TRUE
    )
))

## load data
load('../data/dipperData.RData')

## define hierarchical model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    chi[k] <- 1
    for(j in 1:(k-1))
        chi[k-j] <- (1-phi) + phi * (1-p) * chi[k-j+1]
    for(j in 1:k)
        logChi[j] <- log(chi[j])
    for(i in 1:nind)
        zerosVector[i] ~ dCJS2(length=length[i], nSightingsMinus1=nSightingsMinus1[i], nNonSightings=nNonSightings[i], log_phi=log(phi), log_pp=log(p), log_1minusP=log(1-p), logChi=logChi[last[i]])
})

## define model constants, data, and initial values
last <- apply(y, 1, function(hist) max(which(hist==1)))
nSightings <- apply(y, 1, function(hist) sum(hist, na.rm=TRUE))
constants <- list(k=k, nind=nind, length=last-first, nSightingsMinus1=nSightings-1, nNonSightings=last-first-nSightings+1, last=last)
data <- list(zerosVector=rep(0,nind))
inits <- list(phi=0.6, p=0.9)

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

