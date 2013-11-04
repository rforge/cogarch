\name{COGprm-class}
\Rdversion{1.1}
\docType{class}
\alias{COGprm-class}
\alias{COGprm}
\alias{rCOGARCH-method}

\title{Class \code{"COGprm"}}
\description{
Object class used to define COGARCH parameters. In order to define a generic \eqn{COGARCH(p,q)}, the \code{eta} and \code{phi} can be numeric vectors. With reference to the \eqn{COGARCH(1,1)} definition in Klüppelberg et al. (2004) and Haug et al. (2006), the Slots are defined in order to match the equations
\deqn{dG_t=\sigma_t dL_t}{d G_t=sigma_t d L_t}
\deqn{d\sigma^2_{t"+"}=(\beta-\eta\sigma^2_t)dt+\phi\sigma^2_td\[L,L\]_t^{(d)}}{d sigma^2_{t+}=(beta-eta sigma^2_t)d t+phi sigma^2_t d[L,L]_t^(d)}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{COGprm(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{eta}:}{Object of class \code{"numeric"} auto-regressive parameters.}
    \item{\code{beta}:}{Object of class \code{"numeric"} intercept.}
    \item{\code{phi}:}{Object of class \code{"numeric"} moving-average parameters.}
  }
}
\section{Methods}{
  \describe{
    \item{rCOGARCH}{\code{signature(parameters = "COGprm", increments = "function", obstimes = "numeric", blocks = "numeric", nsim = "numeric", MMM = "numeric", sigmaSq0 = "numeric", L0 = "numeric", Time = "numeric", xname = "character")}: the \code{COGprm} object is used to define the COGARCH parameters where \eqn{L_t} is a Compound Poisson process.}
    \item{rCOGARCH}{\code{signature(parameters = "COGprm", increments = "varGamma", obstimes = "numeric", blocks = "numeric", nsim = "numeric", MMM = "numeric", sigmaSq0 = "numeric", L0 = "numeric", Time = "numeric", xname = "character")}:  the \code{COGprm} object is used to define the COGARCH parameters where \eqn{L_t} is a Variance Gamma process.}
	 }
}
\references{
Klüppelberg, Claudia, Alexander Lindner, and Ross Maller. "A continuous-time GARCH process driven by a Lévy process: Stationarity and second-order behaviour." Journal of Applied Probability 41.3 (2004): 601-622.

Haug, Stephan, Klüppelberg, Claudia, Lindner, Alexander, and Zapp, Mathias. "Method of moment estimation in the COGARCH (1, 1) model." The Econometrics Journal 10.2 (2007): 320-341.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("COGprm")
}
\keyword{classes}
\usage{
COGprm(eta, beta, phi)
}
\arguments{
%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%}
\code{COGprm}-class Slots.
}
