\name{rCOGARCH}
\alias{rCOGARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rCOGARCH(parameters, increments, obstimes, blocks = obstimes, nsim = 1, MMM = 100, sigmaSq0 = 1, L0 = 0, Time = 1, xname = "n", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{increments}{
%%     ~~Describe \code{increments} here~~
}
  \item{obstimes}{
%%     ~~Describe \code{obstimes} here~~
}
  \item{blocks}{
%%     ~~Describe \code{blocks} here~~
}
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
  \item{MMM}{
%%     ~~Describe \code{MMM} here~~
}
  \item{sigmaSq0}{
%%     ~~Describe \code{sigmaSq0} here~~
}
  \item{L0}{
%%     ~~Describe \code{L0} here~~
}
  \item{Time}{
%%     ~~Describe \code{Time} here~~
}
  \item{xname}{
%%     ~~Describe \code{xname} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (parameters, increments, obstimes, blocks = obstimes, 
    nsim = 1, MMM = 100, sigmaSq0 = 1, L0 = 0, Time = 1, xname = "n", 
    ...) 
{
    if (any(length(parameters@phi) > 1, length(parameters@eta) > 
        1)) 
        warning("Only COGARCH(1,1) simulation have been implemented")
    stopifnot(length(parameters@phi) == 1, length(parameters@eta) == 
        1)
    sigmaSq0 <- sigmaSq0
    L0 <- L0
    Time <- Time
    xname <- xname
    MMM <- MMM
    blocks <- blocks
    if (MMM\%\%1 != 0) {
        warning("MMM is rounded to the closest integer")
        MMM <- round(MMM)
    }
    dt <- 1/MMM
    Time <- blocks[1]
    shock <- standardGeneric("rCOGARCH")
    res <- rLatentCOGARCH(shock@timeinc, shock@inc, p = 1, q = 1, 
        eta = parameters@eta, beta = parameters@beta, phi = parameters@phi, 
        sigmaSq0 = sigmaSq0, L0 = L0)
    nblocks <- length(blocks)
    if (nblocks > 1) {
        nit <- length(res@time)
        res@time <- res@time[c(1, nit)]
        res@sigma <- as.matrix(as.matrix(res@sigma)[c(1, nit), 
            ])
        res@G <- as.matrix(as.matrix(res@G)[c(1, nit), ])
        for (i in 2:nblocks) {
            Time <- blocks[i] - blocks[i - 1]
            shock <- standardGeneric("rCOGARCH")
            tmp <- rLatentCOGARCH(shock@timeinc, shock@inc, p = 1, 
                q = 1, eta = parameters@eta, beta = parameters@beta, 
                phi = parameters@phi, sigmaSq0 = tail(res@sigma, 
                  1)^2, L0 = tail(res@G, 1))
            res@time <- c(res@time, tail(tmp@time, 1) + res@time[i])
            res@sigma <- as.matrix(rbind(as.matrix(res@sigma), 
                tail(tmp@sigma, 1)))
            res@G <- as.matrix(rbind(as.matrix(res@G), tail(tmp@G, 
                1)))
        }
    }
    res
  }, generic = structure("rCOGARCH", package = "COGARCH"), package = "COGARCH", group = list(), valueClass = character(0), signature = c("parameters", 
"increments", "obstimes", "blocks", "nsim", "MMM", "sigmaSq0", 
"L0", "Time", "xname"), default = `\001NULL\001`, skeleton = (function (parameters, 
    increments, obstimes, blocks = obstimes, nsim = 1, MMM = 100, 
    sigmaSq0 = 1, L0 = 0, Time = 1, xname = "n", ...) 
stop("invalid call in method dispatch to 'rCOGARCH' (no default method)", 
    domain = NA))(parameters, increments, obstimes, blocks, nsim, 
    MMM, sigmaSq0, L0, Time, xname, ...), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
