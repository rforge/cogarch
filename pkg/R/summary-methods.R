setMethod("summary", signature(object="COGARCH"),
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