\name{shock-class}
\Rdversion{1.1}
\docType{class}
\alias{shock-class}
\alias{shock}
\alias{plot,shock,missing-method}
\alias{shrinker,shock,numeric-method}

\title{Class \code{"shock"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{shock(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{timeinc}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{inc}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Time}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "shock", y = "missing")}: ... }
    \item{shrinker}{\code{signature(shock = "shock", dt = "numeric")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
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
showClass("shock")
}
\keyword{classes}
\usage{
shock(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
