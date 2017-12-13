
# convert a color into its transparent equivalent
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


ltr <- function(xmin, xmax, ymin, ymax) {
  alpha <- (ymax-ymin)/(xmax-xmin)
  beta <- ymin-alpha*xmin
  function(x) {
    alpha*x+beta
  }
}

transformFreq <- ltr(0, .6, 0, .8)

transformH2 <- ltr(0, 1, 0, 1)

transformMu <- ltr(1, 6, .7, 1.5)

plotHeritabilities10k <- function(data, ratesContacts, GEValues, transform, initiatePlot=TRUE, transformRC=function(x) 1/x, 
                                  boxwexFactor=1, 
                                  columns=c(), boxlwd=c(), boxwex=c(), col=c(), border=c(), atOffset=c(), pars=list(),
                                  ...) {
  
  # avoid name conflict with column in data.table
  if(initiatePlot) {
    plot(c(), c(), ...)
    #abline(h=seq(0.1, 1, by=.1), col = "lightgray", lty = "dotted", lwd = par("lwd"))
    abline(v=seq(3, 11, by=2), col = "lightgray", lty = "dotted", lwd = par("lwd"))
    #grid(nx=NA, ny=axTicks(2))
  }
  
  for(rc in ratesContacts) {
    if(nrow(data[rateContact==rc])>1) {
      for(i in 1:length(columns)) {
        boxplot(transform(data[rateContact==rc][, columns[i], with=FALSE]), outline=FALSE,
                notch=FALSE, at=transformRC(rc)+atOffset[i]*boxwexFactor, pch=4, col=col[i], border=border[i],
                boxlwd=boxlwd[i], boxwex=boxwex[i]*boxwexFactor, axes=FALSE, pars=pars, add=TRUE)
      }
    }
  }
}

plotTau10k <- function(data, ratesContacts, GEValues, transform, initiatePlot=TRUE, transformRC=function(x) 1/x,  ...) {
  # avoid name conflict with column in data.table
  if(initiatePlot) {
    plot(c(), c(), ...)
    grid()
  }
  for(rc in ratesContacts) {
    if(nrow(data[rateContact==rc])>1) {
      # plot tau
      boxplot(transform(data[rateContact==rc, meanTauD10k]),
              notch=FALSE, at=transformRC(rc)-.25, pch=4, col='magenta', border='magenta', boxwex=.5, outwex=.20, axes=FALSE, add=TRUE)
      # plot tauR
      boxplot(transform(data[rateContact==rc, meanTauR10k]),
              notch=FALSE, at=transformRC(rc)+.25, pch=4, col='blue', border='blue', boxwex=.5, outwex=.20, axes=FALSE, add=TRUE)  
    }
  }
}

# this uses epidemic$counts to plot the distribution in Infected through time
plotTraitThroughTimeZ <- function(epidemic, tStart=1, tEnd=floor(epidemic$time), 
                                  mean=TRUE, median=FALSE, quantiles=list(c('q2.5%', 'q97.5%')), sdRegion=FALSE,
                                  lwdQuantiles=2, colQuantiles=c('darkgrey'), lwdSDRegion=2, colSDRegion='darkgrey',
                                  alpha=NULL, theta=NULL, G0=NULL, t0=NULL, 
                                  meanSampled=NULL, 
                                  xlim=c(0, tEnd-tStart), ylim=c(0, 5), main="", xlab="", ylab="Distribution of z", lastN=5000, ...) {
  
  plot(x=c(), y=c(), xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)
  
  stats <- as.data.table(t(sapply(tStart:tEnd, function(t) {
    pop <- extractPop(epidemic, sampledOnly=TRUE, activeOnly=FALSE, tMax=t, lastN=lastN)
    
    if(is.null(pop)|nrow(pop)==0) {
      c(mean=NA, median=NA, q2.5=NA, q97.5=NA, sd=NA) 
    } else {
      decomposeTrait(pop, GEVs, copy=FALSE)
      pop[, c(mean=mean(z), median=median(z), q2.5=quantile(z, probs=0.025), q97.5=quantile(z, probs=0.975), sd=sd(z))]
    }
  })))
  
  if(!is.null(quantiles)) {
    for(i in 1:length(quantiles)) {
      segments(x0=(tStart:tEnd)-tStart, x1=(tStart:tEnd)-tStart, 
               y0=stats[tStart:tEnd, q2.5], y1=stats[tStart:tEnd, q97.5], 
               col=colQuantiles[i], lwd=lwdQuantiles) 
    }
  }
  
  if(sdRegion) {
    sd <- stats[tStart:tEnd, sd]
    m <- stats[tStart:tEnd, mean]
    segments(x0=(tStart:tEnd)-tStart, x1=(tStart:tEnd)-tStart, 
             y0=m-sd, y1=m+sd, col=colSDRegion, lwd=lwdSDRegion)  
  }
  
  lines(x=(tStart:tEnd)-tStart, y=stats[tStart:tEnd, mean], col='black')
  
  if(!is.null(theta))
    abline(h=theta, col='darkgreen', lwd=.8, lty=2)
  if(!is.null(meanSampled)) {
    abline(h=meanSampled, col='blue', lwd=.8, lty=2)
  }  
  if(!is.null(G0)) {
    for(s in 1:length(t0)) {
      D0 <- sapply(t0[s]:tEnd, function(t) exp(-alpha*(t-t0[s])))
      D1 <- 1-D0
      mu <- D0*G0[s] + D1*theta
      lines(x=(t0[s]:tEnd)-tStart, y=mu, col='green', lty=2)
    }
  }
}

plotDensities <- function(couples, pop, horiz=FALSE, plotByGenotype=TRUE, main='', ...) {
  cols <- c(2:6, 'brown')
  flip <- function(l) {
    x <- l$x
    l$x <- l$y
    l$y <- x
    l
  }

  dens <- density(couples[, z])
  denspop <- density(pop[, z])
  lty <- 2
  if(horiz) {
    dens <- flip(dens)
    denspop <- flip(denspop)
    lty <- 3
  }
  plot(dens, main=main, lty=lty, ...)
  lines(denspop, lty=1, col='black')
  
  xseq <- seq(0, 6, by=.05)
  dnormx <- dnorm(xseq, couples[, mean(z)], couples[, sd(z)])
  dnormpopx <- dnorm(xseq, pop[, mean(z)], pop[, sd(z)])
  if(horiz) {
    lines(x=dnormx, y=xseq, lty=3, col='darkgrey')
    lines(x=dnormpopx, y=xseq, lty=1, col='darkgrey')
    # bugfix: avoid the strange grey horizontal line by overwriting it with a white one
    abline(h=0, col='white')
  } else {
    lines(x=xseq, y=dnormx, lty=2, col='darkgrey')
    lines(x=xseq, y=dnormpopx, lty=1, col='darkgrey')
  }
  
  if(plotByGenotype) {
    for(i in 1:6) {
      try({
        pi1 <- couples[, sum(gD0==i)/length(gD0)]
        densi1 <- density(couples[gD0==i, z])
        densi1$y <- densi1$y*pi1
        lty <- 2
        if(horiz) {
          densi1 <- flip(densi1)
          lty <- 3
        }
        lines(densi1, col=cols[i], lty=lty)
      })
    }  
  }
}

plotDonorRecipRegression <- function(couples, xlim, pointsFrac=.1, H2, ...) {
  cols <- c(2:6, 'brown')
  couples[, col:=cols[gD0]]
  sampPoints <- sample(1:nrow(couples), size=as.integer(nrow(couples)*pointsFrac), replace=FALSE)
  plot(couples[sampPoints, list(zD0, zR0)], pch=20, col=couples[sampPoints, col], type='p', cex=.12,
       xlim=xlim, ylim=xlim, ...)
  
  aggrdata <- couples[, list(m=mean(zR0), q05=quantile(zR0, probs=.05), q95=quantile(zR0, probs=.95), q25=quantile(zR0, probs=.25), q75=quantile(zR0, probs=.75)), keyby=list(zdi=round(zD0/4, 1)*4)]
  
  points(x=aggrdata[, zdi], y=aggrdata[, m], pch=1, col='black', cex=1.2)
  
  ab <- coef(lm(zR0~zD0, data=couples))
  curve(ab[1]+ab[2]*x, xlim[1]-1, xlim[2]+1, lty=4, col='black', add=T, lwd=1.5) 
  
  meand <- couples[, mean(zD0)]
  meanr <- couples[, mean(zR0)]
  aH2 <- c(meanr-H2*meand, H2)
  curve(aH2[1]+aH2[2]*x, xlim[1]-1, xlim[2]+1, lty=4, col='magenta', add=T, lwd=1.5) 
  
  grid()  
  abline(v=couples[, mean(zD0)], lty=2)
  abline(h=couples[, mean(zR0)], lty=3)
  abline(h=couples[, mean(zD0)], lty=2)
}

calcFreqsPopDEnvGene <- function(epidemic) {
  t10k <- max(which(epidemic$counts[, 'nTips']<=10000))
  pop <- extractPop(epidemic, sampledOnly=TRUE, tMax=t10k) 
  
  couples <- extractDRCouples(epidemic, sampledOnly=TRUE, tMax=t10k)
  
  as.vector(t(cbind(pop[, list(count=length(id)), keyby=list(env, gene)][, list(freqPop=count/nrow(pop))],
                    couples[, list(count=length(id)), keyby=list(envd, gD0)][, list(freqD=count/nrow(couples))],
                    couples[, list(count=length(id)), keyby=list(env, gD0)][, list(freqR=count/nrow(couples))])))
}

calcFreqsPopDEnv <- function(epidemic) {
  t10k <- max(which(epidemic$counts[, 'nTips']<=10000))
  pop <- extractPop(epidemic, sampledOnly=TRUE, tMax=t10k) 
  
  couples <- extractDRCouples(epidemic, sampledOnly=TRUE, tMax=t10k)
  
  as.vector(t(cbind(pop[, list(count=length(id)), keyby=list(env)][, list(freqPop=count/nrow(pop))],
                    couples[, list(count=length(id)), keyby=list(envd)][, list(freqD=count/nrow(couples))],
                    couples[, list(count=length(id)), keyby=list(env)][, list(freqR=count/nrow(couples))])))
}

calcFreqsPopDGene <- function(epidemic) {
  t10k <- max(which(epidemic$counts[, 'nTips']<=10000))
  pop <- extractPop(epidemic, sampledOnly=TRUE, tMax=t10k) 
  
  couples <- extractDRCouples(epidemic, sampledOnly=TRUE, tMax=t10k)
  
  as.vector(t(cbind(pop[, list(count=length(id)), keyby=list(gene)][, list(freqPop=count/nrow(pop))],
                    couples[, list(count=length(id)), keyby=list(gD0)][, list(freqD=count/nrow(couples))])))
}

plotDRRegr <- function(eptree, timeStep, GEVs, doPlot=TRUE, doLegend=TRUE, figPart='', mainText="", kappa_1=0, h2stats, pr, 
                       columns, boxlwd, boxwex, col, border, atOffset, pars, ...) {
  eptreelist <- eptree[, ll.(epidemic, tree10k, list(.[[1]], .[[2]]))]
  couples <- rbindlist(lapply(eptreelist, function(.) {
    if(!is.null(.[[1]]) & !is.null(.[[2]]) & length(.[[2]]$tip.label) > 1000) {
      extractDRCouples(.[[1]], ids <- .[[2]]$tip.label)
    } else {
      NULL
    }
  }))
  
  reportCouples <- b(data=couples, GEValues=GEVs, report=TRUE)
  couples <- reportCouples$couples
  
  nCouples <- nrow(couples)
  
  pop <- rbindlist(lapply(eptreelist, function(.) {
    extractPop(.[[1]], ids <- .[[2]]$tip.label)
  }))
  
  H2aov <- rA(data=pop, GEValues=GEVs, report=TRUE)
  R2adj <- R2adj(data=pop, GEValues=GEVs)
  
  pop <- decomposeTrait(pop, GEVs)
  varCompPop <- decomposeVar(pop)
  
  GValues <- pop[, list(G=unique(G)), keyby=list(gtype=gene)]
  EValues <- pop[, list(E=unique(E)), keyby=list(envtype=env)]
  IValues <- pop[, list(I=unique(I)), keyby=list(envtype=env, gtype=gene)]
  GEValues <- data.table(envtype=c(rep(1, 6), rep(2, 6)), gtype=c(1:6, 1:6), GE=c(GEVs[1, ], GEVs[2, ]), key=c("envtype", "gtype"))

  pop[, eSp:=eSpec(e, gene)]
  
  nPop <- nrow(pop)
  
  couples[, idCouple:=1:nrow(couples)]
  
  idDonor <- couples[, idD]
  idRecip <- couples[, id]

  couples[, ED:=EValues[list(envd), E]]
  couples[, ER:=EValues[list(env), E]]
  
  couples[, GD0:=GValues[list(gD0), G]]
  couples[, ID0:=IValues[list(envd, gD0), I]]
  couples[, eD0estim:=zD0-GD0-ED-ID0]
  
  couples[, GR0:=GValues[list(gD0), G]]
  couples[, IR0:=IValues[list(env, gR0), I]]
  couples[, eR0estim:=zR0-GR0-ER-IR0]
  
  couples[, GD:=GValues[list(gD), G]]
  couples[, ID:=IValues[list(envd, gD), I]]
  couples[, eDestim:=zD-GD-ED-ID]
  
  couples[, GR:=GValues[list(gD), G]]
  couples[, IR:=IValues[list(env, gR), I]]
  couples[, eRestim:=zR-GR-ER-IR]
  
  freqs.ygD0 <- couples[, list(freqs=length(zD0)/nCouples), keyby=list(envd, gD0)][, freqs]
  freqs.ygR0 <- couples[, list(freqs=length(zD0)/nCouples), keyby=list(env, gD0)][, freqs]
  
  freqs.ygp <- pop[, list(freqs=length(z)/nPop), keyby=list(env, gene)][, freqs]
  
  mu <- pop[, mean(z)]
  mu.d <- couples[, mean(zD0)]
  mu.r <- couples[, mean(zR0)]
  
  mu.eSp <- pop[, mean(eSp)]
  sigma2.eSp <- pop[, var(eSp)]
  pVal.eSp <- t.test(pop[, eSp])$p.value
  
  mu.eR0 <- couples[, mean(eR0)]
  sigma2.eR0 <- couples[, var(eR0)]
  pVal.eR0 <- t.test(couples[, eR0])$p.value
  
  mu.eD0 <- couples[, mean(eD0)]
  sigma2.eD0 <- couples[, var(eD0)]
  pVal.eD0 <- t.test(couples[, eD0])$p.value
  
  sigmaDR0 <- couples[, cov(zD0, zR0)]
  sigmaGDGR0 <- couples[, cov(GD0, GR0)]
  sigmaIDIR0 <- couples[, cov(ID0, IR0)]
  sigmaRDRR0 <- couples[, cov(zD0-GD0, zR0-GR0)]
  sigmaeDeR0 <- couples[, cov(eD0, eR0)]
  sigmaGEdGEr <- couples[, cov(GED0, GER0)]
  sigmaeD0eR0estim <- couples[, cov(eD0estim, eR0estim)]
  
  sigmaEDER0 <- couples[, cov(ED, ER)]
  
  sigma2d <- couples[, var(zD0)]
  sigma2r <- couples[, var(zR0)]
  
  stats <- c(mu=mu, mu.d=mu.d, mu.r=mu.r, 
             sigmaDR0=sigmaDR0, 
             sigmaGDGR0=sigmaGDGR0, sigmaEDER0=sigmaEDER0, sigmaIDIR0=sigmaIDIR0, sigmaRDRR0=sigmaRDRR0, sigmaGEdGEr=sigmaGEdGEr, 
             sigmaeD0eR0estim=sigmaeD0eR0estim, sigmaeDeR0=sigmaeDeR0,
             sigma2d=sigma2d, sigma2r=sigma2r, beta0=sigmaDR0/sigma2d,  
             mu.eSp=mu.eSp, sigma2.eSp=sigma2.eSp, pVal.eSp=pVal.eSp, mu.eR0=mu.eR0, sigma2.eR0=sigma2.eR0, pVal.eR0=pVal.eR0,
             mu.eD0=mu.eD0, sigma2.eD0=sigma2.eD0, pVal.eD0=pVal.eD0)
  
  round <- function(x, minDigits) {
    digits <- ifelse(x<0, minDigits, max(minDigits, min(-log10(x), 4)))
    base::round(x, digits)
  }
  
  text <- rep(expression(paste(" ","")), 21)
  text[[1]] <- bquote(paste("",1/kappa==.(kappa_1), ", #all: ", .(nPop), ", #DR: ", .(nCouples), ""))
  text[[2]] <- bquote(paste(" ", {H^2==R[adj]^2}==.(round(R2adj, 2)) ))
  
  #text[[3]] <- bquote(paste(" ", {s^2}(z)==.(round(varCompPop['varz'], 2)), ", ", {s^2}(z[d])==.(round(stats['sigma2d'], 2)), ", ", {s^2}(z[r])==.(round(stats['sigma2r'], 2))))
  text[[3]] <- bquote(paste(" ")) 
  
  text[[4]] <- bquote(paste(" ", {s^2}(G)==.(round(varCompPop['varG'], 2)), ", ", {s^2}(e)==.(round(varCompPop['varz']-varCompPop['varG'], 2)), ", ", {s^2}(z)==.(round(varCompPop['varz'], 2))))
  #text[[4]] <- bquote(paste(" ", {s^2}(G)==.(round(varCompPop['varG'], 2)), ", ", 
  #                          {s^2}(I)==.(round(varCompPop['varI'], 2)), ", ", 
  #                          {s^2}(E)==.(round(varCompPop['varE'], 2)), ", ",
  #                          {s^2}(e)==.(round(varCompPop['varEpsilon'], 2)), ", "))
  # text[[4]] <- bquote(paste(" ", {s^2}(G)==.(round(varCompPop['varG'], 2))))
  
  #text[[5]] <- bquote(paste(" ", -2*{s}(G,E)==.(round(varCompPop['covGE'], 4))))
  text[[5]] <- bquote(paste(" "))
  
  text[[6]] <- bquote(paste(" ", b[0]==.(round(couples[, cov(zD0, zR0)/var(zD0)], 2)), ", ", s(z["d,0"],z["r,0"])==.(round(couples[, cov(zD0, zR0)], 2)), ", ", {s^2}(z["d,0"])==.(round(couples[, var(zD0)], 2))))
  
  
  
  # covariance decomposition at moment of infection (see commented out lines)
  # text[[6]] <- bquote(paste(" ", b[0]==.(round(couples[, cov(zD0, zR0)/var(zD0)], 2)), ", ", s(z["d,0"],z["r,0"])==.(round(couples[, cov(zD0, zR0)], 2)),
  #                           " : ", s(G[d],G[r])==.(round(couples[, cov(GD0, GR0)], 2)), ", "))
  
  # text[[7]] <- bquote(paste("    ",
  #                           s(G[d], I[r]+E[r]+e[r])==(.(round(couples[, cov(GD0, IR0)], 2)))+(.(round(couples[, cov(GD0, ER)], 2)))+(.(round(couples[, cov(GD0, eR0estim)], 2))) ))
  text[[7]] <- bquote(paste(" "))
  # text[[8]] <- bquote(paste("    ",
  #                           s(I[d]+E[d]+e[d], G[r])==(.(round(couples[, cov(ID0, GR0)], 2)))+(.(round(couples[, cov(ED, GR0)], 2)))+(.(round(couples[, cov(eD0estim, GR0)], 2))) ))
  text[[8]] <- bquote(paste(" ", b[tau]==.(round(couples[, cov(zD, zR)/var(zD)], 2)), ", ", s(z[d],z[r])==.(round(couples[, cov(zD, zR)], 2)), ", ", {s^2}(z[d])==.(round(couples[, var(zD)], 2)) ))
  
  # covariance decomposition at moment of sampling (see commented out lines)
  # text[[9]] <- bquote(paste(" ", b[tau]==.(round(couples[, cov(zD, zR)/var(zD)], 2)), ", ", s(z[d],z[r])==.(round(couples[, cov(zD, zR)], 2)),
  #                           " : ", s(G[d],G[r])==.(round(couples[, cov(GD, GR)], 2)), ", "))
  text[[9]] <- bquote(paste("    "))
  
  text[[10]] <- bquote(paste(" ", "ANOVA: ", r[A],"[id]=",.(round(H2aov$H2aov, 2)),
                             ", ", {hat(sigma)^2}(G)==.(round(H2aov$sigmaG2, 2)), 
                             ", ", {hat(sigma)^2}(e)==.(round(H2aov$sigmae2, 2)) ))
  
  # text[[10]] <- bquote(paste("    ",
  #                            s(G[d], I[r]+E[r]+e[r])==(.(round(couples[, cov(GD, IR)], 2)))+(.(round(couples[, cov(GD, ER)], 2)))+(.(round(couples[, cov(GD, eRestim)], 2))) ))
  text[[11]] <- bquote(paste("    "))
  # text[[11]] <- bquote(paste("    ",
  #                            s(I[d]+E[d]+e[d], G[r])==(.(round(couples[, cov(ID, GR)], 2)))+(.(round(couples[, cov(ED, GR)], 2)))+(.(round(couples[, cov(eDestim, GR)], 2))) ))
  text[[12]] <- bquote(paste(" ", "PMM: ", bar(H[BMe]^2)==.(round(h2stats[process==pr&rateContact==1/kappa_1, mean(H2.BMe)], 2)),
                             ", ", bar(sigma[e]^2)==.(round(h2stats[process==pr&rateContact==1/kappa_1, mean(varE.BM)], 2))))
  
  text[[13]] <- bquote(paste("    "))
  
  text[[14]] <- bquote(paste(" ", "POUMM: ", bar(H[OUe]^2)==.(round(h2stats[process==pr&rateContact==1/kappa_1, mean(H2.OUe)], 2)),
                             ", ", bar(sigma[e]^2)==.(round(h2stats[process==pr&rateContact==1/kappa_1, mean(varE.OU)], 2)) ))
  text[[15]] <- bquote(paste(" Mean values: ", bar(z)==.(round(stats['mu'], 2)), ", ", bar(z[d])==.(round(stats['mu.d'], 2)), ", ", bar(z[r])==.(round(stats['mu.r'], 2)))) 
  text[[16]] <- bquote(paste("    ", bar(G[d])==.(round(couples[, mean(GD0)], 2)), ", ", 
                             bar(e[d])==.(round(couples[, mean(eD0estim)], 2)), ", ", 
                             bar(G[r])==.(round(couples[, mean(GR0)], 2)), ", ", 
                             bar(e[r])==.(round(couples[, mean(eR0estim)], 2))))
  text[[17]] <- bquote(paste(" GE-values: y=1:      ", .(toString(round(GEVs[1, 1:6], 2)))))
  text[[18]] <- bquote(paste("                    y=2:        ", .(toString(round(GEVs[2, 1:6], 2)))))
  text[[19]] <- " Environment type X genotype frequencies:"
  text[[20]] <- bquote(paste(" population: y=1: ", .(toString(round(sum(freqs.ygp[1:6]), 2))), " | ", .(toString(round(freqs.ygp[1:6], 2)))))
  text[[21]] <- bquote(paste("                   y=2: ", .(toString(round(sum(freqs.ygp[7:12]), 2))), " | ", .(toString(round(freqs.ygp[7:12], 2)))))
  text[[22]] <- bquote(paste(" donors:      y=1: ", .(toString(round(sum(freqs.ygD0[1:6]), 2))), " | ", .(toString(round(freqs.ygD0[1:6], 2)))))
  text[[23]] <- bquote(paste("                   y=2: ", .(toString(round(sum(freqs.ygD0[7:12]), 2))), " | ", .(toString(round(freqs.ygD0[7:12], 2)))))
  text[[24]] <- bquote(paste(" recipients:  y=1: ", .(toString(round(sum(freqs.ygR0[1:6]), 2))), " | ", .(toString(round(freqs.ygR0[1:6], 2)))))
  text[[25]] <- bquote(paste("                   y=2: ", .(toString(round(sum(freqs.ygR0[7:12]), 2))), " | ", .(toString(round(freqs.ygR0[7:12], 2)))))
  
  if(doPlot) {
    factor=2
    par(cex=1/factor)
    par(cex.axis=0.8*factor)
    par(cex.lab=0.9*factor)
    par(cex.main=1*factor)
    par(tcl=-0.5)
    par(mgp=c(1.6, .5, 0)*factor)
    # avoid bold titles
    par(font.main=1)
    
    par(mar=c(0, 5, 6, 0))
    
    couples[, z:=zD0]
    couples[, envi:=envd]
    
    plotDensities(couples, pop, plotByGenotype=FALSE, 
                  main=mainText,
                  xlim=c(-0.2,  6.2), ylim=c(0, .8), xaxt='n', yaxt='n', 
                  ylab='f', xlab='', frame=FALSE)
    
    
    mtext(text=text[[1]], side=3, line=0.3, cex=.7, adj=0, at=-0.5)
    axis(side=2, at=seq(0, .8, by=.2))
    
    transf <- ltr(0, 0.5, 0.2, 0.8)
    plotHeritabilities10k(data=h2stats[process==pr], ratesContacts=1/kappa_1,
                          GEValues=GEVs, 
                          transform=transf, initiatePlot=FALSE, 
                          transformRC=function(rc) 5.8, boxwexFactor=.7, 
                          columns=columns, boxlwd=boxlwd, boxwex=boxwex, col=col, 
                          border=border, atOffset=atOffset, pars=pars)
    labs=seq(0, .5, by=.1)
    axis(side=2, pos=5.1, at=transf(labs), labels=labs)
    mtext(text=expression(H^2), line=-5.2, at=4.1, cex=.8)
    
    mtext(text=figPart, font=2, side=3, line=2.5, cex=1.1, adj=0, at=-2.1)
    
    if(doLegend) {
      legend(x='topleft', legend=c(expression(f(z)), expression(f(z[d])), 
                                   expression(paste("N[ ", bar(z), ", ", {s^2}(z), "]")), 
                                   expression(paste("N[ ", bar(z[d]), ", ", {s^2}(z[d]), "]"))), 
             lty=c(1, 2, 1, 2), col=c(1, 1, 'darkgrey', 'darkgrey'), cex=1, bg='white')
    }
    par(mar=c(5, 5, 0, 0))
    plotDonorRecipRegression(couples, xlim=c(-0.2,  6.2), xlab=expression(z[d]), ylab=expression(z[r]), 
                             H2=varCompPop['varG']/varCompPop['varz'], ...)
    if(doLegend) {
      legend(x='topleft', legend=c(expression(bar(z[d])), 
                                   expression(bar(z[r])),
                                   expression(hat(z)["r,OLS"](z[d])), 
                                   expression({{hat(z)[paste("r,", H^2)]}}(z[d])),
                                   expression({bar(z[r])}[paste(" | ",z[d])])), 
             lty=c(2, 3, 4, 4, NA), col=c(1, 1, 1, 'magenta', 1), pch=c(NA, NA, NA, NA, 1), cex=1, bg='white')
    }
    par(mar=c(0, 0.4, 6, 0))
    
    plot(c(), frame=FALSE, xlim=c(1, 13), ylim=c(-.3, .5), xlab='', ylab='', xaxt='n', yaxt='n')
    
    
    line <- rbind(c(1, 2,   3,   4,      5, 6,    7,   8,    9,    10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,   22,   23,   24, 25),
                  c(1, 2.25, 3.5, 4.5, 5.5, 7, 8.15, 9.3, 10.8, 11.95, 13.1, 14.4,   16, 17.6, 18.6, 19.7, 20.7, 21.7, 22.7, 23.7, 24.7, 25.7, 26.7, 27.7, 28.7))
    col <- c(rep('black', 16), 'black', 'darkgray', 'black', rep(c('black', 'darkgrey'), 3))
    cex <- c(0.7, 0.6, rep(.5, 23))
    mtext(text=text[2:14], line=2.2-line[2,2:14], cex=cex[2:14],
          col=col[2:14], adj=0)
  
    couples[, z:=zR0]
    couples[, envi:=env]
    
    par(mar=c(5, 0, 0, 6))
    plotDensities(couples, pop, plotByGenotype=FALSE, horiz=TRUE, xlim=c(0, .8), ylim=c(-0.2,  6.2), xaxt='n', yaxt='n', ylab='', xlab='f', frame=FALSE)
    axis(side=1, at=seq(0, .8, by=.2))
    #grid()
    if(doLegend) {
      legend(x='bottomright', legend=c(expression(f(z)), expression(f(z[r])), 
                                       expression(paste("N[ ", bar(z), ", ", {s^2}(z), "]")), 
                                       expression(paste("N[ ", bar(z[r]), ", ", {s^2}(z[r]), "]"))), 
             lty=c(1, 3, 1, 3), col=c('black', 'black', 'darkgrey', 'darkgrey'), cex=1, bg='white')
    }
    
  }
  list(pop=pop, couples=couples, varCompPop=varCompPop, stats=stats, GValues=GValues, EValues=EValues, IValues=IValues, GEValues=GEValues)
}

plotMeans <- function(data, ratesContacts, GEValues, transform, ...) {
  # avoid name conflict with column in data.table
  GEVs <- GEValues
  for(rc in ratesContacts) {
    if(nrow(data[rateContact==rc])>1) {
      # plot empirical mus
      boxplot(transform(data[rateContact==rc, ll.(mlOUTree, tree, nTips, 
                                              expr=if(is.list(.[[1]]))
                                                calcMuAtTime(attr(.[[1]]$value, 'grmax'), 
                                                             .[[1]]$par[1], .[[1]]$par[2], .[[1]]$par[3], quantile(nodeTimes(.[[2]])[1:.[[3]]], probs=1)/100) 
                                              else NA)]),
              at=1/rc-.12, col='green', border='green', boxwex=.3, axes=F, add=T)
      boxplot(transform(sim[, mean]), at=1/rc+.12, col='black', border='black', boxwex=.3, axes=F, add=T)
    }
  }
}
