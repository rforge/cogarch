setMethod("plot",
          signature(x = "shock", y = "missing"),
          function (x, y, ...) 
          {
            matplot(cumsum(x@timeinc), apply(x@inc, 2, cumsum), type="s", xlab="t", ylab="", xlim=c(0, x@Time), ...)
          }
)