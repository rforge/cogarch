\name{COGprm}
\alias{COGprm}
\alias{COGprm-class}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Class For Cogarch Parameters}
\description{
The \code{COGprm} class is a class of the \code{COGARCH} package that contains the values of the Cogarch(1,1) model.

The Cogarch Model is defined as:

\code{dG_{t} = sigma_{t-}dL_{t}}

\code{dsigma_{t}^{2}=(beta-eta*sigma^{2}_{t-})dt+phi*sigma^{2}_{t-}d[L]^{d}_{t}}

where \code{beta > 0}, \code{eta > 0}  and \code{eta > 0}.

\code{[L]^{d}_{t}} denotes the discrete part of the quadratic variation of the driving Levy process \code{L_{t}
}.
}
\usage{
COGprm(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{Arguments for the constructor \code{COGprm(eta,beta,phi)}:}
  
  \describe{
  \item{eta}{is a numeric object.}
  \item{beta}{is a numeric object.}
  \item{phi}{is a numeric object.}
  }
}


\details{
An object of class \code{COGprm} is the main input for the \code{\link{rCOGARCH}} method 
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Klupperlberg, C., Lindner, A., and Maller, R (2004). \emph{A continuous-time GARCH process driven by a Levy process: stationarity and second-order behaviour}. J. Appl. Probab., 41(3) 601-622.
}

\author{
Cogarch Team
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
