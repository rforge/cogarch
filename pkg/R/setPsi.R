setGeneric("setPsi",function(cogpar,levy){
	standardGeneric("setPsi")
	})
	
setMethod("setPsi", c(cogpar="COGprm", levy="varGamma"),
          function(cogpar, levy) {
          	beta<-cogpar@beta
			eta<-cogpar@eta
			phi<-cogpar@phi
			C<-levy@c
			tau<-0
            list(psi1= -eta +phi, psi2= -2*eta+2*phi+3*phi^2/C, psi3= -3*(eta-phi-10*phi^3/C^2-3*phi^2/C), psi4= -4*eta + (2*phi*(2*C^3 + 9*C^2*phi + 60*C*phi^2 + 315*phi^3))/C^3, tau=tau)
          })
