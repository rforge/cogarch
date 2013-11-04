setGeneric("shrinker",
           function(shock, dt)
             standardGeneric("shrinker")
)

setMethod("shrinker", c("shock", "numeric"),
          function(shock, dt) {
            m <- 1/dt
            timeinc <- shock@timeinc
            timeinc[length(timeinc)+1] <- shock@Time - sum(timeinc)
            repetitions <- round(m*timeinc)
            repetitions[repetitions==0] <- 1
            deltat <- rep(timeinc/repetitions, repetitions)
            tot <- length(deltat)
            inc <- as.matrix(rep(0, tot))
            idx <- cumsum(repetitions)
            idx <- idx[-length(idx)]
            inc[idx, 1] <- shock@inc[, 1]
            new("shock", timeinc=deltat, inc=inc, Time=sum(deltat))
          }
)