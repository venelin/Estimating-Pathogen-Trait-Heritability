# Here we setup the toy-model, specifying the pathogen strains (here, called genotypes), 
# the general host-type x carried-strain effects (GEVs), 
# the SIR parameters (birth/death/recovery/transmission/mutation rate), 
# the total population sizes and the time for simulation

# genotype encodings with their allele contents
genotypes <- rbind(c(0, 0, 1, 0, 1),
                   c(0, 0, 1, 1, 0),
                   c(0, 1, 0, 0, 1),
                   c(0, 1, 0, 1, 0),
                   c(1, 0, 0, 0, 1),
                   c(1, 0, 0, 1, 0))

# indices for genotypes
genotypes <- cbind(1:6, genotypes)
colnames(genotypes) <- c('gindex', 'x13', 'x12', 'x11', 'x22', 'x21')

# data.table provides a nice syntax here
genotypes <- as.data.table(genotypes)

# add the contents of allelic combinations
genotypes[, xx1322:=x13*x22]
genotypes[, xx1321:=x13*x21]
genotypes[, xx1222:=x12*x22]
genotypes[, xx1221:=x12*x21]
genotypes[, xx1122:=x11*x22]
genotypes[, xx1121:=x11*x21]


# convert back to a matrix
genotypes <- as.matrix(genotypes)

# combinations of genotypes and environments
genotypeXenv <- cbind(gi=rep(genotypes[, 1], 2), 
                      ei=c(rep(1, 6), rep(2, 6)), 
                      gxei=1:12,
                      rbind(genotypes[, -1], genotypes[, -1]), 
                      y2=c(rep(0, 6), rep(1, 6)), 
                      y1=c(rep(1, 6), rep(0, 6)))


# make results reproducible (using same seed as for the discrete generation case)
set.seed(5)

pe <- c(.5, .5)
n <- 2

# General host-type x strain effects (expected phenotypes for genotype by environment combinations)
GEVs <- normReactPlot(genotypeXenv, runif(n=12, min=2, max=4)[c(7:12, 1:6)], do.plot=TRUE)

sigmae <- .6

timeStep <- 0.05

pg.init <- c(1, 0, 0, 0, 0, 0)
sde <- c(sigmae, sigmae)

# death-rate as a function of viral load and natural death rate mu
rateDie <- function(z, mu) {
  V <- 10^z
  Dmin <- 2
  Dmax <- 25*12
  D50 <- 10^3
  Dk <- 1.4
  (V^Dk+D50^Dk)/(Dmin*(V^Dk+D50^Dk)+((Dmax-Dmin)*D50^Dk)) + mu
}

meanRateDie <- function(z, mu) {
  rep(0.01, length(z))
}

# infectionrate as a function of viral load and rate of risky contacts
rateInfect <- function(z, rateContact) {
  V <- 10^z
  Emin <- .3
  Emax <- .6
  E50 <- 10^3
  Ek <- 1.4
  E <- Emin+(Emax-Emin)*V^Ek/(V^Ek+E50^Ek)
  E*rateContact
}

meanRateInfect <- function(z, rateContact) {
  rep(0.45*rateContact, length(z))
}

# per locus mutation rates
rateMutate <- function(GEValues, es, envs, genes) {
  z <- es+GEValues[cbind(envs, genes)]
  V <- 10^z
  Mmin <- 0.00
  Mmax <- 0.2
  M50 <- 10^3
  Mk <- 1.4
  Mmin+(Mmax-Mmin)*V^Mk/(V^Mk+M50^Mk)
}

meanRateMutate <- function(GEValues, es, envs, genes) {
  rep(0.01, length(genes))
}

rateSampleNeutral <- 1/(4*12)

# maximum time before starting graceful fadeout of the epidemic
maxTime <- 2400

# continue the epidemic outbreak until reaching ...
maxNTips <- 10000

# is the special environmental effect unique for each pathogen genoetype in an individual
eUniqForEachG <- TRUE

# are only beneficial (i.e. increasing the trait-value) mutations allowed
selectWithinHost <- TRUE

# time to continue the simulation of transmission after reaching maxNTips (this was introduced
# in order to study post-outbreak dynamics, i.e. epidemic waves after exhaustion of the 
# susceptible pool)
expandTimeAfterMaxNTips <- 4

#initial population size (equilibrium in the absence of disease)
N <- 1e5

# setting the between-host and within-host dynamics of the simulation:

# natural death rate
mu <- 1/850

# constant birth rate that maintains this equilibrium
nu <- ifelse(is.finite(N), mu*N, 0)

# multiplier for within-host transition rates (see function simulateEpidemic in the patherit package)
rateTransTemplate_32 <- function() {
  rbind(
    c(1, 1/2, 0, 1/2, 0),
    c(1, 0, 1/2, 0, 1/2),
    c(1/2, 0, 1, 1/2, 0),
    c(0, 1/2, 1, 0, 1/2),
    c(1/2, 0, 1/2, 0, 1),
    c(0, 1/2, 0, 1/2, 1))
}

# Simulate one epidemic and fit POUMM to the population of the first 10000 cases
simulateAndFit10k <- function(p, doTreeAndFit=TRUE) {
  set.seed(p$id*29)
  if(is.null(p$epidemic)) {
    epidemic <- do.call(simulateEpidemic, list(Ninit=p$N, nu=p$nu, mu=p$mu, 
                                               pe=p$pe, sde=p$sde, pg.init=p$pg.init,
                                               GEValues=p$GEValues, rateContact=p$rateContact, rateInfect=p$rateInfect, 
                                               rateDie=p$rateDie, rateSample=p$rateSample, 
                                               rateMutate=p$rateMutate, rateTransTemplate=p$rateTransTemplate,
                                               eUniqForEachG=p$eUniqForEachG, selectWithinHost=p$selectWithinHost, timeStep=p$timeStep, 
                                               maxTime=p$maxTime, maxNTips=p$maxNTips, 
                                               expandTimeAfterMaxNTips=p$expandTimeAfterMaxNTips, process=p$process))
    if(doTreeAndFit) {
      treeAll <- extractTree(epidemic)
    } else {
      treeAll <- NULL
    }
  } else {
    epidemic <- p$epidemic
    treeAll <- p$treeAll
  }
  
  if(doTreeAndFit & !is.null(treeAll)) {
    vAll <- epidemic$gen[treeAll$tip.label, calcValue(env, gene, e, p$GEValues)]
    tipTimesAll <- nodeTimes(treeAll)[1:length(treeAll$tip.label)]
    
    timeThr10k <- sort(tipTimesAll)[10000]
    if(!is.na(timeThr10k)) {
      tree10k <- extractTree(epidemic, tips=which(tipTimesAll<=timeThr10k)) 
      v10k <- vAll[which(tipTimesAll<=timeThr10k)]
      nTips10k <- length(tree10k$tip.label)  
    } else {
      tree10k <- treeAll
      v10k <- vAll
      nTips10k <- length(tree10k$tip.label)
    }
    
    # fit POUMM and PMM using the default patherit settings.
    anH2 <- estimateH2(v10k, tree10k, 
                       methods = list(PP=TRUE, 
                                      POUMM=list(nSamplesMCMC = 0, verbose=TRUE), 
                                      PMM=list(nSamplesMCMC = 0, verbose=TRUE)),   
                       verbose = TRUE, usempfr=-2)
    
    corrProfile <- corrProfile(anH2)
    # # code using the old patherit version 
    # mlOUTree10k <- ml.poumm(v10k, tree10k, distgr='maxlik', parMax=c(alpha=31, theta=10, sigma=20, sigmae=10), divideEdgesBy=100)
    # mlBMTree10k <- ml.poumm(v10k, tree10k, distgr='maxlik', parMax=c(sigma=20, sigmae=10), parFixed=c(alpha=0, theta=0), divideEdgesBy=100)
  } else {
    vAll <- NULL
    tree10k <- NULL
    v10k <- NULL
    nTips10k <- 0
    # mlOUTree10k <- NULL
    # mlBMTree10k <- NULL
  }
  
  list(epidemic=epidemic, treeAll=treeAll, vAll=vAll, tree10k=tree10k, v10k=v10k, 
       nTips10k=nTips10k, 
       #mlOUTree10k=mlOUTree10k, mlBMTree10k=mlBMTree10k,
       anH2 = anH2, 
       corrProfile = corrProfile
       )
}