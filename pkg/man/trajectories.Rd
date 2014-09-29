\name{trajectories}
\alias{trajectories}
\alias{trajectories-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Class for sample paths of a Cogarch(1,1) process}
\description{An object of class \code{trajectories} contains the Cogarch(1,1) sample paths.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trajectories(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{\dots}{Arguments for the constructor \code{trajectories(time,sigma,G)}:
\describe{
  \item{time}{is a \code{numeric} object containing the observation dates of the process.}
  \item{sigma}{is a \code{matrix} object in which the sample paths of the volatility process have been stored.}
  \item{G}{is a \code{matrix} object containing the simulated trajectories of the Cogarch(1,1).}
  }
  }
}
\details{The class is built internally by the function \code{\link{rCOGARCH}} and contains the simulated sample path for a Cogarch(1,1) model}
\author{
Cogarch Team
}
\section{Methods}{ 
\describe{ 
  \item{plot}{Plot methods} 
  \item{summary}{Summary methods}
  } 
}