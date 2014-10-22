setMethod("as.vector", signature(x="COGprm"),
          function(x) {theta<-c(x@beta,x@eta,x@phi)
          names(theta)<-c("beta","eta","phi")
          theta}
          )