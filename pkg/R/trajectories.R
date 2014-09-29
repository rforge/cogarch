trajectories <- setClass(
  "trajectories",
  representation(
    time = "numeric",
    sigma = "matrix",
    G = "matrix"
  )
)

# Plot

setMethod("plot",
          signature(x = "trajectories", y = "missing"),
          function (x, y, ...) 
          {
            par(mfrow=c(2,1))
            matplot(x@time, x@G, type="s", xlab="t", ylab="", xlim=c(0, max(x@time)), main="G", ...)
            matplot(x@time, x@sigma, type="s", xlab="t", ylab="", xlim=c(0, max(x@time)), main="sigma", ...)
          }
)

# Summary
setMethod("summary", signature(object="trajectories"),
          function(object) {
            minT <- min(object@time)
            maxT <- max(object@time)
            avglag <- mean(diff(object@time))
            nobs <- nrow(object@G)
            nsim <- ncol(object@G)
            cat("\nSIMULATIONS OF COGARCH(1,1) PROCESS\n",
                "\nnumber of trajectories:  ", nsim,
                "\nnumber of observations:  ", nobs,
                "\ninitial time:            ", minT,
                "\nterminal time:           ", maxT,
                "\naverage observations lag:", avglag)
          }
)
