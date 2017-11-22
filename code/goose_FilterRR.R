## Filter (RR) MCMC for the goose model

## load nimble library
library(nimble)

## define custom distribution
dDHMM <- nimbleFunction(
    run = function(x = double(1), length = double(), prior = double(1), condition = double(1), Z = double(3), useZt = double(), T = double(3), useTt = double(), mult = double(), log = double()) {
        pi <- prior
        logL <- 0
        Zind <- 1
        Tind <- 1
        for(t in 1:length) {
            if(useZt) Zind <- t
            Zcurrent <- Z[,,Zind]
            if(t == 1) {
                for(i in 1:dim(Zcurrent)[1])
                    Zcurrent[i, ] <- Zcurrent[i, ] * condition[i]
                for(j in 1:dim(Zcurrent)[2]) {
                    s <- sum(Zcurrent[ ,j])
                    if(s != 0) Zcurrent[ ,j] <- Zcurrent[ ,j] / s
                }
            }
            Zpi <- Zcurrent[x[t], ] * pi
            sumZpi <- sum(Zpi)
            logL <- logL + log(sumZpi) * mult
            if(t != length) {
                if(useTt) Tind <- t
                pi <- (T[,,Tind] %*% asCol(Zpi) / sumZpi)[ ,1]
            }
        }
        returnType(double())
        if(log) return(logL) else return(exp(logL))
    }
)

rDHMM <- nimbleFunction(
    run = function(n = integer(), length = double(), prior = double(1), condition = double(1), Z = double(3), useZt = double(), T = double(3), useTt = double(), mult = double()) {
        declare(x, double(1, length))
        returnType(double(1))
        return(x)
    }
)

registerDistributions(list(
    dDHMM = list(
        BUGSdist = 'dDHMM(length, prior, condition, Z, useZt, T, useTt, mult)',
        types = c('value = double(1)', 'length = double()', 'prior = double(1)', 'condition = double(1)', 'Z = double(3)', 'useZt = double()', 'T = double(3)', 'useTt = double()', 'mult = double()'),
        discrete = TRUE
    )
))

## load data
load('../data/gooseData.RData')

## define hierarchical model
code <- nimbleCode({
    for(i in 1:6)
        p[i] ~ dunif(0, 1)
    for(t in 1:3) {
        Z[1,1,t] <- p[1]
        Z[2,1,t] <- 0
        Z[3,1,t] <- 0
        Z[4,1,t] <- 1 - p[1]
        Z[1,2,t] <- 0
        Z[2,2,t] <- p[2]
        Z[3,2,t] <- 0
        Z[4,2,t] <- 1 - p[2]
        Z[1,3,t] <- 0
        Z[2,3,t] <- 0
        Z[3,3,t] <- p[3]
        Z[4,3,t] <- 1 - p[3]
        Z[1,4,t] <- 0
        Z[2,4,t] <- 0
        Z[3,4,t] <- 0
        Z[4,4,t] <- 1
    }
    Z[1,1,4] <- p[4]
    Z[2,1,4] <- 0
    Z[3,1,4] <- 0
    Z[4,1,4] <- 1 - p[4]
    Z[1,2,4] <- 0
    Z[2,2,4] <- p[5]
    Z[3,2,4] <- 0
    Z[4,2,4] <- 1 - p[5]
    Z[1,3,4] <- 0
    Z[2,3,4] <- 0
    Z[3,3,4] <- p[6]
    Z[4,3,4] <- 1 - p[6]
    Z[1,4,4] <- 0
    Z[2,4,4] <- 0
    Z[3,4,4] <- 0
    Z[4,4,4] <- 1
    for(i in 1:3)
        phi[i] ~ dunif(0, 1)
    for(i in 1:2)
        for(j in 1:3)
            for(l in 1:2) {
                alpha[i,j,l] ~ dnorm(0, sd = 100)
                exal[i,j,l] <- exp(alpha[i,j,l])
            }
    for(j in 1:3)
        for(l in 1:2) {
            psi[1,j,l] <- exal[1,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[2,j,l] <- exal[2,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[3,j,l] <-           1 / (1 + exal[1,j,l] + exal[2,j,l])
        }
    for(t in 1:2)
        for(j in 1:3) {
            for(i in 1:3)
                T[i,j,t] <- phi[j] * psi[i,j,1]
            T[4,j,t] <- 1 - phi[j]
        }
    for(t in 3:4)
        for(j in 1:3) {
            for(i in 1:3)
                T[i,j,t] <- phi[j] * psi[i,j,2]
            T[4,j,t] <- 1 - phi[j]
        }
    for(t in 1:4) {
        T[1,4,t] <- 0
        T[2,4,t] <- 0
        T[3,4,t] <- 0
        T[4,4,t] <- 1
    }
    for (i in 1:nind) {
        y[i, first[i]:k] ~ dDHMM(length = k-first[i]+1, prior = prior[1:4], condition = condition[1:4], Z = Z[1:k,1:k,first[i]:k], useZt = 1, T = T[1:k,1:k,first[i]:k], useTt = 1, mult = mult[i])
    }
})

## define model constants, data, and initial values
constants <- list(nind=nind, k=k, first=first, mult=mult, prior=c(1/3, 1/3, 1/3, 0), condition = c(1, 1, 1, 0))
data <- list(y=y)
inits <- list(p=rep(1/2,6), phi=rep(1/2,3), alpha=array(0,c(2,3,2)))

## create model object
Rmodel <- nimbleModel(code, constants, data, inits, check = FALSE)

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

