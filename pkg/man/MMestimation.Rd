\name{MMestimation}
\alias{MMestimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameter estimation by Method of Moments (MM)
%%  ~~function to do ... ~~
}
\description{
Parameter estimation by Method of Moments (MM) as described in the reference below
}
\usage{
MMestimation(x, d = sqrt(length(x)), explicit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
is a vector of equispaced log returns (according to the model the are increments of the G process)
}
  \item{d}{
an integer that defines how many observations in the past are used 
}
  \item{explicit}{
logical (default is \code{TRUE}), if false minimization is performed numerically
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The value is a list with two entries, \code{theta} that is a COGprm object containing the estimated parameters and \code{aux} that contains auxiliary informations

}
\references{
Haug, S., Kl{\"u}ppelberg, C., Lindner, A., and Zapp, M. (2007).
\emph{Method of moment estimation in the} COGARCH(1,1) \emph{model.}  Econom. J., 10(2):320--341.
}
\author{
Cogarch Team
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(COGARCH)
test.COGprm.object<-COGprm(eta=0.053, beta=0.04, phi=0.038)
test.VG.object<-varGamma(c=1)
set.seed(123)
traj<-rCOGARCH(test.COGprm.object, test.VG.object,1:10000, nsim=2, dt=0.01)
test.trajectory<-diff(traj@G[,1])
MMestimation(test.trajectory, d=70, explicit = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
