\name{rLatentCOGARCH}
\alias{rLatentCOGARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rLatentCOGARCH(deltat, inc, p = 1, q = 1, eta, beta, phi, sigmaSq0, L0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deltat}{
%%     ~~Describe \code{deltat} here~~
}
  \item{inc}{
%%     ~~Describe \code{inc} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{sigmaSq0}{
%%     ~~Describe \code{sigmaSq0} here~~
}
  \item{L0}{
%%     ~~Describe \code{L0} here~~
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
function (deltat, inc, p = 1, q = 1, eta, beta, phi, sigmaSq0, 
    L0 = 0) 
{
    n <- length(deltat)
    stopifnot(sigmaSq0 > 0)
    B <- beta * deltat
    A <- 1 - eta * deltat + phi * inc^2
    X_t <- apply(A, 2, cumprod)
    X_s <- apply(beta * deltat/X_t, 2, cumsum)
    sigmaSq <- rbind(sigmaSq0, X_t * (X_s + matrix(sigmaSq0, 
        nrow(X_s), ncol(X_s), byrow = TRUE)))
    sigma <- sqrt(sigmaSq)
    infraret <- rbind(as.matrix(inc), 0) * sigma
    G <- apply(rbind(L0, as.matrix(infraret[-n, ])), 2, cumsum)
    res <- new("COGARCH", time = c(0, cumsum(deltat)), sigma = sigma, 
        G = G)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
