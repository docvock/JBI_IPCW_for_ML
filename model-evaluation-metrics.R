###################################
# File: ./model-evaluation-metrics.R
#		Purpose: Create functions which return the cNRI (censored NRI), Harrell's C-index, and calibration
#							test statistic and make the calibration plots.
#   Author: David M. Vock
#   Last Modified: January 8, 2016
#   Relies On: none
#   Files Created: none
####################################


## Function for computing net reclassification improvement, taking censoring into account
## Due to Pencina (2011)
## Typically cNRI weights reclassification for events and non-events equally but this need not 
## be the case. This is especially problematic when the proportion of events is significantly
## different from 0.5 and one of the models is misspecified. The argument class_wt gives the
## relative weight for reclassifying events and non-event respectively

cNRI <- function(p1, p2, cutpts, test.dat, t, class_wt = c(1, 1)) {
  risks.1 <- cut(1-p1,breaks=cutpts,labels=FALSE)
  risks.2 <- cut(1-p2,breaks=cutpts,labels=FALSE)
  T <- test.dat$T.use
  C <- test.dat$C.use
  
  n <- length(T)

  up.class <- (risks.2>risks.1)
  down.class <- (risks.2<risks.1)
  n.u <- sum(up.class) ## Number up-classified
  n.d <- sum(down.class) ## Number down-classified

  KM <- survfit(Surv(T,C)~1) ## overall
  p.KM <- 1 - KM$surv[max(which(KM$time<=t))] ## P(event)
  
  if(n.u==0) {
    p.KM.u <- 1
    } else {
      KM.u <- survfit(Surv(T,C)~1,subset=up.class) ## up-classified
      p.KM.u <- 1 - KM.u$surv[max(which(KM.u$time<=t))] ## P(event|up)
    }
  if(n.d==0) {
    p.KM.d <- 1
  } else {
    KM.d <- survfit(Surv(T,C)~1,subset=down.class) ## down-classified
    p.KM.d <- 1 - KM.d$surv[max(which(KM.d$time<=t))] ## P(event|down)
  }

  nri.e <- (n.u*p.KM.u - n.d*p.KM.d)/(n*p.KM)
  nri.ne <- (n.d*(1-p.KM.d) - n.u*(1-p.KM.u))/(n*(1-p.KM))
	wt.e <- class_wt[1]
	wt.ne <- class_wt[2]
  
  c(cNRI.events=nri.e,cNRI.nonevents=nri.ne,cNRI=wt.e*nri.e+wt.ne*nri.ne)
}

## Harrell's C-index
Cindex <- function(p, test.dat) {
  test.dat$p <- 1-p
  Cindex_all <- survConcordance(Surv(T.use, C.use) ~ p, data = test.dat)
  Cindex <- Cindex_all$concordance
  return(Cindex)
}

## Calibration statistic 
calib.stat <- function(p, cutpts, test.dat, verbose = FALSE) {
	risk.class <- cut(p, cutpts, labels=FALSE)
  lev.stats <- sapply(1 : (length(cutpts) - 1), function(f) {
  	ind <- which(risk.class == f)
  	if(length(ind)>0){
  	S.KM <- survfit(Surv(test.dat$T.use[ind], test.dat$C.use[ind]) ~ 1, error = "greenwood")
  	ind.surv <- max(which(S.KM$time <= (SURVTIME*365)))
  	p.KM <- S.KM$surv[ind.surv]
  	S.KM.std.err <- S.KM$std.err[ind.surv]*p.KM
    if (verbose == TRUE) {
    	print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE), S.KM.std.err,
    		(mean(p[ind], na.rm = TRUE) - p.KM)^2 / (S.KM.std.err)^2))
    }
    (mean(p[ind], na.rm = TRUE) - p.KM)^2 / (S.KM.std.err)^2
  }else {
  	0
  }})
  calib.stat <- (sum(lev.stats, na.rm=TRUE))
  return(calib.stat)
}

## Calibration plot 
calib.plot <- function(p.total, cutpts, test.dat, main.title, pch.use = 19, col.use = "black",
	lty.use = 1, legend.use = FALSE, legend.text = NULL, legend.loc = "bottomleft") {
	if (length(pch.use == 1)) {
		pch.use <- rep(pch.use, length(p.total))
	}
	if (length(col.use == 1)) {
		col.use <- rep(col.use, length(p.total))
	}
	if (length(lty.use == 1)) {
		lty.use <- rep(lty.use, length(p.total))
	}
	t1 <- lapply(p.total, calib.func, cutpts=cutpts, test.dat=test.dat)
	#calib.data <- calib.func(p, cutpts, test.dat)
	resid.tot <- NULL
	obs.tot <- NULL
	for (j in 1:length(p.total)) {
		resid.tot <- c(resid.tot, t1[[j]][,1] - t1[[j]][,2])
		obs.tot <- c(obs.tot, t1[[j]][,2])	
	}
	plot(t1[[1]][,2], t1[[j]][,1] - t1[[j]][,2], type = "n", xlim = c(min(obs.tot), max(obs.tot)), 
		ylim = c(min(resid.tot), max(resid.tot)), xlab = "Observed Risk", 
		ylab = "Predicted - Observed Risk", main = main.title )
	for (j in 1:length(p.total)) {
		resid <- t1[[j]][,1] - t1[[j]][,2]
		obs <-  t1[[j]][,2]
		lines(obs, resid, type = "b", pch = pch.use[j], col = col.use[j], lty = lty.use[j], lwd = 3)
	}
	abline(h = 0, col = "red", lwd = 3)
	if (legend.use == TRUE) {
		legend(legend.loc, legend.text, pch = pch.use, col = col.use, lty = lty.use)
	}
}

calib.func <- function(p, cutpts, test.dat) {
	risk.class <- cut(p, cutpts, labels=FALSE)
	lev.stats <- sapply(1 : (length(cutpts) - 1), function(f) {
		ind <- which(risk.class == f)
		S.KM <- calcSurv(Surv(test.dat$T.use[ind], test.dat$C.use[ind]))
		ind.surv <- max(which(S.KM$t <= (SURVTIME*365)))
		p.KM <- S.KM$SKM[ind.surv]
		c(mean(p[ind], na.rm = TRUE),p.KM)})
	
	calib.func <- 1-t(lev.stats)
	colnames(calib.func) <- c("Predicted", "Observed")
	return(calib.func)
}
