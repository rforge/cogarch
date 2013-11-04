rComPoisson <- function(Time, expr = rep(1, n), lambda=1, xname = "n") {  
  
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    expr <- call(as.character(sexpr), as.name(xname))
  }
  else {
    if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in% 
            all.vars(sexpr))) 
      stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'", 
                    xname), domain = NA)
    expr <- sexpr
  }
  
  n <- 0
  deltat <- NULL
  while(sum(deltat) <= Time) {
    n <- n + 1
    deltat[n] <- rexp(1, 1/lambda)
  }
  deltat <- deltat[-n]
  n <- n - 1 # numero di salti, nell'ultimo intervallo non ci sono salti ma serve per definire
  # quando termina il processo
  ll <- list(n = n)
  names(ll) <- xname
  j <- eval(expr, envir = ll, enclos = parent.frame())
  
  res <- new("shock", timeinc = deltat, inc = as.matrix(j), Time=Time)
  return(res)
}