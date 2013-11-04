setGeneric("rCOGARCH",
           function(
             parameters,
             increments,
             obstimes,
             blocks = obstimes,
             nsim = 1,
             MMM = 100,
             sigmaSq0 = 1,
             L0 = 0,
             Time = 1,
             xname = "n",
             ...) {
             
             if (any(length(parameters@phi)>1, length(parameters@eta)>1)) warning("Only COGARCH(1,1) simulation have been implemented")
             stopifnot(length(parameters@phi)==1, length(parameters@eta)==1)
             
             sigmaSq0 <- sigmaSq0
             L0 <- L0
             Time <- Time
             xname <- xname
             MMM <- MMM
             blocks <- blocks
             
             if (MMM %% 1 != 0) {
               warning("MMM is rounded to the closest integer")
               MMM <- round(MMM)
             }
             dt <- 1/MMM
             
             Time <- blocks[1]
             shock <- standardGeneric("rCOGARCH")
             res <- rLatentCOGARCH(shock@timeinc, shock@inc, p=1, q=1, eta=parameters@eta, beta=parameters@beta, phi=parameters@phi, sigmaSq0=sigmaSq0, L0=L0)
             
             nblocks <- length(blocks)
             if (nblocks > 1) {
               nit <- length(res@time)
               res@time <- res@time[c(1, nit)]
               res@sigma <- as.matrix(as.matrix(res@sigma)[c(1, nit),])
               res@G <- as.matrix(as.matrix(res@G)[c(1, nit),])
               for(i in 2:nblocks) {
                 Time <- blocks[i] - blocks[i-1]
                 shock <- standardGeneric("rCOGARCH")
                 tmp <- rLatentCOGARCH(shock@timeinc, shock@inc, p=1, q=1, eta=parameters@eta, beta=parameters@beta, phi=parameters@phi, sigmaSq0=tail(res@sigma, 1)^2, L0=tail(res@G, 1))
                 res@time <- c(res@time, tail(tmp@time, 1) + res@time[i]) # ultimo elemento o ultima riga
                 res@sigma <- as.matrix(rbind(as.matrix(res@sigma), tail(tmp@sigma, 1)))
                 res@G <- as.matrix(rbind(as.matrix(res@G), tail(tmp@G, 1)))
               }
             }
             res
           }
)

setMethod("rCOGARCH", c("COGprm", "varGamma", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character"),
          function(parameters, increments, obstimes, blocks, nsim, MMM, sigmaSq0, L0, Time, xname) {
            shape <- increments@delta*increments@c
            rate <- sqrt(2*increments@c)
            incshock <- sapply(1:nsim, function(i) rgamma(Time*MMM, shape, rate) - rgamma(Time*MMM, shape, rate))
            new(Class = "shock", timeinc=rep(1/MMM, Time*MMM), inc=as.matrix(incshock), Time=Time)
          })

setMethod("rCOGARCH", c("COGprm", "function", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character"),
          function(parameters, increments, obstimes, blocks, nsim, MMM, sigmaSq0, L0, lambda, Time, xname) { # non mettiamo il default a lambda cosi' l'utente deve ricordarsi di passarlo
            
            dt <- 1/MMM
            stopifnot(is.numeric(lambda))
            shock <- lapply(1:nsim, function(i) {
              shrinker(rComPoisson(Time=Time, expr=increments, lambda=lambda, xname=xname), dt)
            })
            timeshock <- sort(unique(unlist(lapply(shock, function(x) cumsum(x@timeinc)))))
            funshock <- lapply(shock, function(x) approxfun(c(0, cumsum(x@timeinc), Time), c(0, cumsum(as.vector(x@inc)), sum(as.vector(x@inc))), method="constant")) # ATTENZIONE! lavoro a una traiettoria alla volta
            incshock <- as.data.frame(lapply(funshock, function(x) x(timeshock)))
            incshock <- as.matrix(incshock)
            incshock <- diff(incshock)
            colnames(incshock) <- paste("T", 1:ncol(incshock), sep="")
            new("shock", timeinc=diff(timeshock), inc=incshock, Time=Time)
          })