varGamma <- setClass("varGamma", representation(
  #  delta="numeric", # e' il dt... lo controlliamo noi o lasciamo che l'utente possa sbagliare?
  c = "numeric"
))

setGeneric("rCOGARCH",
           function(
             parameters,
             increments,
             obstimes,
             blocks = obstimes,
             nsim = 1,
             dt = 1/100,
             sigmaSq0 = 1,
             G0 = 0,
             Time = 1,
             xname = "n",
             ...) {
               
             MMM <- round(1/dt)
             
             if (any(length(parameters@phi)>1, length(parameters@eta)>1)) warning("Only COGARCH(1,1) simulation has been implemented")
             stopifnot(length(parameters@phi)==1, length(parameters@eta)==1)
             
             sigmaSq0 <- sigmaSq0
             G0 <- G0
             # Time <- Time
             xname <- xname
             dt <- dt
             nsim <- nsim
             
             blocks <- blocks
             
             Time <- blocks[1]
             shock <- standardGeneric("rCOGARCH")
             res <- rLatentCOGARCH(shock@timeinc, shock@inc, p=1, q=1, eta=parameters@eta, beta=parameters@beta, phi=parameters@phi, sigmaSq0=sigmaSq0, G0=G0)
             
             nblocks <- length(blocks)
             if (nblocks > 1) {
               dates <- res@time
               wch <- c(which(dates %in% obstimes), length(res@time))
               wch <- unique(wch)
               res@time <- res@time[c(1, wch)]
               res@sigma <- as.matrix(as.matrix(res@sigma)[c(1, wch),])
               res@G <- as.matrix(as.matrix(res@G)[c(1, wch),])
               for(i in 2:nblocks) {
                 Time <- blocks[i] - blocks[i-1]
                 shock <- standardGeneric("rCOGARCH")
                 tmp <- rLatentCOGARCH(shock@timeinc, shock@inc, p=1, q=1, eta=parameters@eta, beta=parameters@beta, phi=parameters@phi, sigmaSq0=tail(res@sigma, 1)^2, G0=tail(res@G, 1))
#                  wch <- ((tmp@time + res@time[i]) %in% obstimes)[-1]
                 dates <- tmp@time + tail(res@time, 1)
                 wch <- which(dates %in% obstimes)
                 wch <- unique(wch)[-1]
                 wch <- unique(c(wch, length(tmp@time)))
                 res@time <- c(res@time, dates[wch]) # ultimo elemento o ultima riga
                # res@sigma <- as.matrix(rbind(as.matrix(res@sigma), as.matrix(as.matrix(tmp@sigma)[wch,])))
                res@sigma <- as.matrix(rbind(as.matrix(res@sigma), t(as.matrix(as.matrix(tmp@sigma)[wch,]))))

                 # res@G <- as.matrix(rbind(as.matrix(res@G), as.matrix(as.matrix(tmp@G)[wch,])))
                
                res@G <- as.matrix(rbind(as.matrix(res@G), t(as.matrix(as.matrix(tmp@G)[wch,]))))
               }
             }
             #wch <- which(res@time %in% obstimes) #
             wch<-approx(x=obstimes,y=obstimes, xout=res@time,  method="const")  
             wch2<- which(unique(wch$y) %in% obstimes)
             res@time <- res@time[na.omit(wch2)]
             res@time <- c(0,res@time)
             res@sigma <- as.matrix(as.matrix(res@sigma)[na.omit(wch2),])
             res@sigma <- rbind(matrix(sigmaSq0,nrow=1,ncol=nsim),res@sigma)
             res@G <- as.matrix(as.matrix(res@G)[na.omit(wch2),])
             res@G <- rbind(matrix(G0,nrow=1,ncol=nsim),res@G)
             res
           }
)

setMethod("rCOGARCH", c("COGprm", "varGamma", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character"),
          function(parameters, increments, obstimes, blocks, nsim, dt, sigmaSq0, G0, Time, xname) {
            shape <- dt*increments@c
            rate <- sqrt(2*increments@c)
            incshock <- sapply(1:nsim, function(i) rgamma(Time/dt, shape, rate) - rgamma(Time/dt, shape, rate))
            new(Class = "shock", timeinc=rep(dt, Time/dt), inc=as.matrix(incshock), Time=Time)
          })

setMethod("rCOGARCH", c("COGprm", "function", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character"),
          function(parameters, increments, obstimes, blocks, nsim, dt, sigmaSq0, G0, lambda, Time, xname) { # non mettiamo il default a lambda cosi' l'utente deve ricordarsi di passarlo
            
            MMM <- round(1/dt)
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