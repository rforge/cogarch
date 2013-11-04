\name{COGARCH-package}
\alias{COGARCH-package}
\docType{package}
\title{
COGARCH simulation
}
\description{
The package provides simulation S4 procedure for the COGARCH introduced by Klüppelberg et al. (2004)
}
\details{
\tabular{ll}{
Package: \tab COGARCH\cr
Type: \tab Package\cr
Version: \tab 0.1 \cr
Date: \tab 2013-10-03\cr
License: \tab GPL (>= 2) \cr
Depends: \tab methods\cr
}
This package contains S4 functions that provide an efficient simulation method for COGARCH based on the solution of stochastic recurrence equations.
}
\author{
A, Arcagni
L, Mercuri

Maintainer: A, Arcagni <alberto.arcagni@unimib.it>
}
\references{
Klüppelberg, Claudia, Alexander Lindner, and Ross Maller. "A continuous-time GARCH process driven by a Lévy process: Stationarity and second-order behaviour." Journal of Applied Probability 41.3 (2004): 601-622.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{package; COGARCH; simulation; estimation}

\examples{
# number of simulated values per time unit
MMM <- 100
# COGARCH(1, 1) parameters
prm <- COGprm(eta = 0.053, beta = 0.04, phi = 0.038)
# parameters of the Variance Gamma increment for the underlying Lévy process
prmVgamma <- varGamma(delta = 1/MMM, c = 1)
# simulation of 10 trajectories
prv <- rCOGARCH(prm, prmVgamma, obstimes = 1:2000/10, nsim = 10, sigmaSq0 = 0.2, L0 = 0, MMM = MMM)
plot(prv)
}
