# RUNNING BISSE-NESS --------------------------------------------

## Due to a change in sample() behaviour in newer R it is necessary to
## use an older algorithm to replicate the previous examples
if (getRversion() >= "3.6.0") {
  RNGkind(sample.kind = "Rounding")
}

# forcing states to be binary
st <- setNames(alldata$tridepth, alldata$sp)
st[st == "shallow"] <- 0
st[st %in% c("deep", "intermediate")] <- 1
mode(st) <- "numeric"


## This builds the likelihood of the data according to BiSSEness:
tree <- force.ultrametric(tree)
lik <- make.bisseness(tree, st)


## ML search:  First we make hueristic guess at a starting point, based
## on the constant-rate birth-death model assuming anagenesis (uses
## \link{make.bd})
startp <- starting.point.bisse(tree)


## We then take the total amount of anagenetic change expected across
## the tree and assign half of this change to anagenesis and half to
## cladogenetic change at the nodes as a heuristic starting point:
t <- branching.times(tree)
tryq <- 1/2 * startp[["q01"]] * sum(t)/length(t)
p <- c(startp[1:4], startp[5:6]/2, p0c=tryq, p0a=0.5, p1c=tryq, p1a=0.5)



## Start an ML search from this point.  
fit <- find.mle(lik, p, method="subplex")

## Compare the fit to a constrained model that only allows the trait
## to change along a lineage (anagenesis).  
lik.no.clado <- constrain(lik, p0c ~ 0, p1c ~ 0)
fit.no.clado <- find.mle(lik.no.clado,p[argnames(lik.no.clado)])


cat("Results from ML estimates of anagenesis vs. cladogenesis: ")
anova(fit, no.clado=fit.no.clado)



# this takes a LONG time
if(rerun){
# MCMC run: We use the ML estimate from the full model
## as a starting point.
ml.start.pt <- pmax(coef(fit), 1e-4)
# have to tweak numbers because running into rounding issues
ml.start.pt[ml.start.pt > 0.9] <- 0.89

## Make exponential priors for the rate parameters and uniform priors
## for the cladogenetic change probability prarameters.
make.prior.exp_ness <- function(r, min=0, max=1) {
  function(pars) {
    sum(dexp(pars[1:6], rate=r, log=TRUE)) +
      sum(dunif(pars[7:10], min, max, log=TRUE))
  }
}


## Choosing the slice sampling parameter, w (affects speed):
library(numDeriv)
hess <- hessian(lik, ml.start.pt)
vcv <- -solve(hess)
sehess <- sqrt(abs(diag(vcv)))
w <- 2 * pmin(sehess, .2)

## Setting the priors
r <- log(length(tree$tip.label))/max(branching.times(tree))
prior <- make.prior.exp_ness(1/(2*r))
prior(ml.start.pt)

## Running the mcmc chain
steps <- 1000
nchains <- 1
for (i in 1:nchains){
  output <- mcmc(lik, ml.start.pt, nsteps=steps, w=w, prior=prior) 
  saveRDS(output, here(paste0("data/prog_data/bisseness_output_", suffix, "_", nchains, ".Rdata")))
}
}


### plotting MCMCs
output <- readRDS(here(paste0("data/prog_data/bisseness_output_", suffix, "_", nchains, ".Rdata")))
## Plotting one mcmc run (dropping first 500 points)
mcmcdrop <- subset(output, i > 100)
mcmcdrop$r0 <- mcmcdrop$lambda0 - mcmcdrop$mu0
mcmcdrop$r1 <- mcmcdrop$lambda1 - mcmcdrop$mu1



## Plotting the effect of the state on diversification and character; code from Magnuson-Ford et al. 2012 supplemental materials

## change parameters
col.fill <- c("#FF000066", "#0000FF66")
col.line <- c("red", "blue")
ylab <- "Probability density"
op <- par(mfrow = c(1,2), mar=c(4.1, 4.5, .5, .5), oma=c(0, 0, 1.5, 0))
# op <- par(mfcol=c(3,2), mar=c(4.1, 4.5, .5, .5), oma=c(0, 0, 1.5, 0))
# profiles.plot(mcmcdrop[c("lambda0", "lambda1")], col.line, col.fill,
#               las=1, xlab="Speciation rate", ylab=ylab)
# legend("topright", c("Photic (State 0)", "Aphotic (State 1)"),
#        fill=col.fill)
# curve(dexp(x, 1/(2 * r)), add=TRUE)
# 
# profiles.plot(mcmcdrop[c("mu0", "mu1")], col.line, col.fill,
#               las=1, xlab="Extinction rate", ylab=ylab)
# curve(dexp(x, 1/(2 * r)), add=TRUE)
# 
# profiles.plot(mcmcdrop[c("r0", "r1")], col.line, col.fill,
#               las=1, xlab="Diversification rate", ylab=ylab)

profiles.plot(mcmcdrop[c("q01", "q10")], col.line, col.fill,
              las=1, xlab="Anagenetic transition rate", ylab=ylab)
legend("topright", c("Photic", "Aphotic"), fill=col.fill, bty = "n")
curve(dexp(x, 1/(2 * r)), add=TRUE)

profiles.plot(mcmcdrop[c("p0c", "p1c")], col.line, col.fill,
              las=1, xlab="Cladogenetic transition rate", ylab="")
curve(dunif(x), add=TRUE)

# profiles.plot(mcmcdrop[c("p0a", "p1a")], col.line, col.fill,
#               las=1, xlab="Prob. assymetric change", ylab=ylab)
# curve(dunif(x), add=TRUE)
