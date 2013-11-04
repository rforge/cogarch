\name{rComPoisson}
\alias{rComPoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rComPoisson(Time, expr = rep(1, n), lambda = 1, xname = "n")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Time}{
%%     ~~Describe \code{Time} here~~
}
  \item{expr}{
%%     ~~Describe \code{expr} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{xname}{
%%     ~~Describe \code{xname} here~~
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
function (Time, expr = rep(1, n), lambda = 1, xname = "n") 
{
    sexpr <- substitute(expr)
    if (is.name(sexpr)) {
        expr <- call(as.character(sexpr), as.name(xname))
    }
    else {
        if (!((is.call(sexpr) || is.expression(sexpr)) && xname \%in\% 
            all.vars(sexpr))) 
            stop(gettextf("'expr' must be a function, or a call or an expression containing '\%s'", 
                xname), domain = NA)
        expr <- sexpr
    }
    n <- 0
    deltat <- NULL
    while (sum(deltat) <= Time) {
        n <- n + 1
        deltat[n] <- rexp(1, 1/lambda)
    }
    deltat <- deltat[-n]
    n <- n - 1
    ll <- list(n = n)
    names(ll) <- xname
    j <- eval(expr, envir = ll, enclos = parent.frame())
    res <- new("shock", timeinc = deltat, inc = as.matrix(j), 
        Time = Time)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
