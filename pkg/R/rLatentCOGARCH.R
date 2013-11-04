rLatentCOGARCH <- function(deltat, inc, p=1, q=1, eta, beta, phi, sigmaSq0, L0=0) {
  n <- length(deltat)
  stopifnot(sigmaSq0 > 0)
  B <- beta*deltat
  A <- 1 - eta*deltat + phi*inc^2
  
  X_t <- apply(A, 2, cumprod) # processo ausiliario
  X_s <- apply(beta*deltat/X_t, 2, cumsum)
  sigmaSq <- rbind(sigmaSq0, X_t*(X_s+matrix(sigmaSq0, nrow(X_s), ncol(X_s), byrow=TRUE)))
  sigma <- sqrt(sigmaSq)
  infraret <- rbind(as.matrix(inc), 0)*sigma
  
  G<-apply(rbind(L0, as.matrix(infraret[-n,])), 2, cumsum)
  
  res <- new("COGARCH", time=c(0, cumsum(deltat)), sigma=sigma, G=G)
}