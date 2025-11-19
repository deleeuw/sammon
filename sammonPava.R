sammonPava <- function(y) {
  nobj <- length(y)
  amat <- matrix(0, nobj, nobj)
  diag(amat) <- 1
  amat[outer(1:nobj, 1:nobj, function(i, j) (i - j) == 1)] <- -1
  bmat <- solve(amat)
  xold <- sqrt(y)
  mumu <- sum(xold)
  xold <- xold / mumu
  sold <- sum(y / xold)
  repeat {
    lbd <- bmat %*% 
  }
}