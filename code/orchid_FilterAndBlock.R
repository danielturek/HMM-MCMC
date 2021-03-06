## Filter and Block MCMC for the orchid model

## load nimble library
library(nimble)

## define custom distribution
dDHMMorchid <- nimbleFunction(
    run = function(x = double(1), length = double(), prior = double(1), Z = double(2), T = double(3), log = double()) {
        pi <- prior
        logL <- 0
        for(t in 1:length) {
            Zpi <- Z[x[t], ] * pi
            sumZpi <- sum(Zpi)
            logL <- logL + log(sumZpi)
            if(t != length)   pi <- (T[,,t] %*% asCol(Zpi) / sumZpi)[ ,1]
        }
        returnType(double())
        return(logL)
    }
)

rDHMMorchid <- nimbleFunction(
    run = function(n = integer(), length = double(), prior = double(1), Z = double(2), T = double(3)) {
        declare(x, double(1, length))
        returnType(double(1))
        return(x)
    }
)

registerDistributions(list(
    dDHMMorchid = list(
        BUGSdist = 'dDHMMorchid(length, prior, Z, T)',
        types = c('value = double(1)', 'length = double()', 'prior = double(1)', 'Z = double(2)', 'T = double(3)'),
        discrete = TRUE,
        mixedSizes = TRUE
    )
))

## load data
load('../data/orchidData.RData')

## define hierarchical model
code <- nimbleCode({
    for (t in 1:(k-1)) {
        s[t] ~ dunif(0, 1)
    }
    for (i in 1:3) {
        a[i] ~ dgamma(1, 1) 
        psiD[i] <- a[i]/sum(a[1:3]) 
        b[i] ~ dgamma(1, 1) 
        psiV[i] <- b[i]/sum(b[1:3]) 
        c[i] ~ dgamma(1, 1) 
        psiF[i] <- c[i]/sum(c[1:3]) 
    }
    for (t in 1:(k-1)) {
        T[1,1,t] <- s[t] * psiV[1]
        T[2,1,t] <- s[t] * psiV[2]
        T[3,1,t] <- s[t] * psiV[3]
        T[4,1,t] <- 1-s[t]
        T[1,2,t] <- s[t] * psiF[1]
        T[2,2,t] <- s[t] * psiF[2]
        T[3,2,t] <- s[t] * psiF[3]
        T[4,2,t] <- 1-s[t]
        T[1,3,t] <- s[t] * psiD[1]
        T[2,3,t] <- s[t] * psiD[2]
        T[3,3,t] <- s[t] * psiD[3]
        T[4,3,t] <- 1-s[t]
        T[1,4,t] <- 0
        T[2,4,t] <- 0
        T[3,4,t] <- 0
        T[4,4,t] <- 1
    }
    T[1,1,k] <- 1
    T[2,1,k] <- 0
    T[3,1,k] <- 0
    T[4,1,k] <- 0
    T[1,2,k] <- 0
    T[2,2,k] <- 1
    T[3,2,k] <- 0
    T[4,2,k] <- 0
    T[1,3,k] <- 0
    T[2,3,k] <- 0
    T[3,3,k] <- 1
    T[4,3,k] <- 0
    T[1,4,k] <- 0
    T[2,4,k] <- 0
    T[3,4,k] <- 0
    T[4,4,k] <- 1
    for (i in 1:nind) {
        y[i, f[i]:k] ~ dDHMMorchid(length = k-f[i]+1, prior = prior[1:4], Z = Z[1:3,1:4], T = T[1:4,1:4,f[i]:k])
    }
})

## define model constants, data, and initial values
constants <- list(f=f, k=k, nind=nind, prior=c(1/2,1/2,0,0), Z=array(c(1,0,0,0,1,0,0,0,1,0,0,1), c(3,4)))
data <- list(y = y)
inits <- list(s = rep(1/2,k-1), a = rep(1,3), b = rep(1,3), c = rep(1,3))

## create model object
Rmodel <- nimbleModel(code, constants, data, inits, check = FALSE)

## specify MCMC algorithm
spec <- configureMCMC(Rmodel, autoBlock = TRUE)
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

