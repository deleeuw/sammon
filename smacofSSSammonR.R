source("smacofDataUtilities.R")
source("smacofTorgerson.R")
source("smacofPlots.R")

suppressPackageStartupMessages(library(isotone, quietly = TRUE))

delta <- matrix(0, 6, 6)
delta[, 1] <- c(0, 1, 2, 3, 4, 5)
delta[, 2] <- c(0, 0, 6, 7, 8, 9)
delta[, 3] <- c(0, 0, 0, 10, 11, 12)
delta[, 4] <- c(0, 0, 0, 0, 13, 14)
delta[, 5] <- c(0, 0, 0, 0, 0, 15)
delta <- as.dist(delta + t(delta))

smallData <- makeMDSData(delta, weights = NULL)

smacofSSSammonR <- function(theData,
                      ndim = 2,
                      xinit = NULL,
                      ties = 1,
                      itmax = 1000,
                      eps = 1e-6,
                      verbose = FALSE,
                      ordinal = FALSE,
                      weighted = FALSE) {
  
  ties <- switch(ties, "primary", "secondary", "tertiary")
  ndat <- theData$ndat
  nobj <- theData$nobj
  iind <- theData$iind
  jind <- theData$jind
  dhat <- theData$delta
  wght <- theData$weights
  if (!weighted) {
    wght <- rep(1, ndat)
  }
  wght <- wght / sum(wght)
  dhat <- dhat / sum(wght * dhat)
  if (is.null(xinit)) {
    xinit <- smacofTorgerson(theData, ndim)$conf
  }
  xold <- xinit
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    dold[k] <- sqrt(sum((xold[iind[k], ] - xold[jind[k], ])^2))
  }
  labd <- sum(wght * dold) / sum((wght / dhat) * dold^2)
  xold <- labd * xold
  dold <- labd * dold
  sold <- sum(wght * dhat) + sum((wght / dhat) * dold^2) - 2 * sum(wght * dold) 
  itel <- 1
  repeat {
    vmat <- matrix(0, nobj, nobj)
    bmat <- matrix(0, nobj, nobj)
    for (k in 1:ndat) {
      vmat[iind[k], jind[k]] <- -wght[k] / dhat[k]
      bmat[iind[k], jind[k]] <- -wght[k] / dold[k]
    }
    vmat <- vmat + t(vmat)
    bmat <- bmat + t(bmat)
    diag(vmat) <- -rowSums(vmat)
    diag(bmat) <- -rowSums(bmat)
    vinv <- solve(vmat + 1 / nobj) - 1 / nobj
    xmid <- vinv %*% bmat %*% xold
    dmid <- rep(0, ndat)
    for (k in 1:ndat) {
      dmid[k] <- sqrt(sum((xmid[iind[k], ] - xmid[jind[k], ])^2))
    }
    rmid <- dmid / dhat
    smid <- sum(wght * dhat) + sum((wght / dhat) * dmid^2) - 2 * sum(wght * dmid) 
    if (ordinal) {
      dhat <- sqrt(gpava(theData$delta, dmid ^ 2, weights = wght, ties = ties)$x)
      dhat <- dhat / sum(wght * dhat)
      snew <- sum(wght * dhat) + sum((wght / dhat) * dmid^2) - 2 * sum(wght * dmid) 
    } else {
      snew <- smid
    }
    if (verbose) {
      if (ordinal) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "smid ",
          formatC(
            smid,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      } else {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      }
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
    xold <- xmid
    dold <- dmid
  }
  result <- list(
    delta = theData$delta,
    dhat = dhat,
    confdist = dmid,
    conf = xmid,
    weightmat = wght,
    stress = snew,
    ndim = ndim,
    niter = itel,
    nobj = nobj,
    iind = iind,
    jind = jind,
    ordinal = ordinal
  )
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}