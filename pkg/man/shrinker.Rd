\name{shrinker}
\alias{shrinker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shrinker(shock, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shock}{
%%     ~~Describe \code{shock} here~~
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
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
structure(function (shock, dt) 
standardGeneric("shrinker"), generic = structure("shrinker", package = "COGARCH"), package = "COGARCH", group = list(), valueClass = character(0), signature = c("shock", 
"dt"), default = `\001NULL\001`, skeleton = (function (shock, 
    dt) 
stop("invalid call in method dispatch to 'shrinker' (no default method)", 
    domain = NA))(shock, dt), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
