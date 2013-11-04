setMethod("plot",
          signature(x = "COGARCH", y = "missing"),
          function (x, y, ...) 
          {
            par(mfrow=c(2,1))
            matplot(x@time, x@G, type="s", xlab="t", ylab="", xlim=c(0, max(x@time)), main="G", ...)
            matplot(x@time, x@sigma, type="s", xlab="t", ylab="", xlim=c(0, max(x@time)), main="sigma", ...)
          }
)
