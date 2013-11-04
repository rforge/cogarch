\name{COGARCH-class}
\Rdversion{1.1}
\docType{class}
\alias{COGARCH-class}
\alias{COGARCH}
\alias{plot,COGARCH,missing-method}
\alias{summary,COGARCH-method}

\title{Class \code{"COGARCH"}}
\description{
This class is used to contain multiple trajectories of a COGARCH process. Even if the process is continuous, its trajectories are observed at fixed dates reported in the slot \code{time}. The COGARCH trajectories are collected by columns in matrix \code{G}, and the volatility trajectories in matrix \code{sigma}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{COGARCH(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{time}:}{Object of class \code{"numeric"} that contains the fixed dates.}
    \item{\code{sigma}:}{Object of class \code{"matrix"} that contains the volatility trajectories.}
    \item{\code{G}:}{Object of class \code{"matrix"} that contains the COGARCH trajectories.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "COGARCH", y = "missing")}: the method provides a \code{\link{matplot}} of the COGARCH trajectories and a \code{\link{matplot}} of the volatility trajectories setting the slot \code{time} as abscissa.}
    \item{summary}{\code{signature(object = "COGARCH")}: the method provides a short \code{\link{summary} of the trajectories.}
	 }
}
}

\seealso{
\code{\link{rCOGARCH}}, \code{\link{matplot}}, \code{\link{summary}}
}
\examples{
showClass("COGARCH")
}
\keyword{classes}
\usage{
COGARCH(...)
}

\arguments{
  \item{\dots}{
Slots of the class \code{\link{COGARCH-class}}
}
}
\details{
It is not necessary to use the class generator because the trajectories are simulated by the \code{\link{rCOGARCH}} function. To observe the generated trajectories and to have a their description, methods for this class of the functions \code{plot} and \code{summary} are realized.
}
