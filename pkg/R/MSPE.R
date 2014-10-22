setGeneric("MSPEcontrast", function(theta, samp, delta, q.back, ...){
standardGeneric("MSPEcontrast")
}
)


setMethod("MSPEcontrast",c("COGprm", "vector", "numeric", "numeric"), function(theta, samp, delta, q.back,...){#... dovrebbe per ora contenere il vg o altra specifica del levy 
	if (theta@beta>0.000001 && theta@eta>0.000001 &&theta@phi>0.000001){
		if(q.back<length(samp)){
			#delta is related to the fact that observations are from G^(delta)_t=G_t-G_(t-delta) with G_t solving dG_t=sigma_t dL_t
			#I will in part mimic notations in Sorensen's paper
			beta<-theta@beta
			eta<-theta@eta
			phi<-theta@phi
			lista<-setPsi(theta,...)
			psi1<-lista$psi1
			psi2<-lista$psi2
			psi3<-lista$psi3
			psi4<-lista$psi4
			tau<-lista$tau
			samp2<-samp^2   #in input vogliamo i log-returns, non i quadrati, ma se si hanno i log-prezzi bisogna già applicare le diff prima
			lung<-length(samp)
			if (psi4>=0) return(NA)
			apsi1<-abs(psi1)
			apsi2<-abs(psi2)
			EG2<-beta*delta/apsi1
			m<-beta/apsi1
			v<-2*beta^2/abs(psi1*psi2)
			EG4<-6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)+(2*beta^2/phi^2)*(2/apsi2-1/apsi1)*delta+3*m^2*delta^2
			varG2add1<- 6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)
			varG2add2<- 2*delta*beta^2*(2/apsi2-1/apsi1)/phi^2
			varG2add3<-2*m^2*delta^2
			varG2<-varG2add1+varG2add2+varG2add3 # gamma_zero_teorico
			covG2null<- m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(1-exp(-delta*apsi1))*(exp(delta*apsi1)-1)/apsi1
			covG2<-function(h) covG2null*exp(-h*apsi1)
			CC<-matrix(numeric(q.back*q.back),q.back)
			
			CC<-covG2null*exp(-abs(row(CC)-col(CC))*delta*apsi1)
			
			diag(CC)<-numeric(q.back)+varG2
			b<-covG2null*exp(-(1:q.back)*delta*apsi1)
			a<-drop(solve(CC)%*%b)
			a0<-EG2*(1-sum(a))
			calcola.predittore <- function(x,sam) a0+a%*%sam[(x-1):(x-q.back)]
			#confronto.predittore<- function(x,sam) (sam[x]-(a0+a%*%sam[(x-1):(x-q.back)]))^2
			predittori<-sapply((q.back+1):length(samp2),calcola.predittore,sam=samp2)
			#c(sum((samp2[(q.back+1):length(samp2)]-predittori)^2),predittori)
			sum((samp2[(q.back+1):length(samp2)]-predittori)^2)
		}
	else stop("q.back argument needs to be smaller than the samp size")
	}
	else NA
}
)

setMethod("MSPEcontrast",c("vector", "vector", "numeric", "numeric"), function(theta, samp, delta, q.back,...){#... dovrebbe per ora contenere il vg o altra specifica del levy 
	if (theta[1]>0.000001 && theta[2]>0.000001 &&theta[3]>0.000001){
		if(q.back<length(samp)){
			#delta is related to the fact that observations are from G^(delta)_t=G_t-G_(t-delta) with G_t solving dG_t=sigma_t dL_t
			#I will in part mimic notations in Sorensen's paper
			beta<-theta[1]
			eta<-theta[2]
			phi<-theta[3]
			lista<-setPsi(COGprm(eta=eta, beta=beta, phi=phi),...)
			psi1<-lista$psi1
			psi2<-lista$psi2
			psi3<-lista$psi3
			psi4<-lista$psi4
			tau<-lista$tau
			samp2<-samp^2   #in input vogliamo i log-returns, non i quadrati, ma se si hanno i log-prezzi bisogna già applicare le diff prima
			lung<-length(samp)
			if (psi4>=0) return(NA)
			apsi1<-abs(psi1)
			apsi2<-abs(psi2)
			EG2<-beta*delta/apsi1
			m<-beta/apsi1
			v<-2*beta^2/abs(psi1*psi2)
			EG4<-6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)+(2*beta^2/phi^2)*(2/apsi2-1/apsi1)*delta+3*m^2*delta^2
			varG2add1<- 6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)
			varG2add2<- 2*delta*beta^2*(2/apsi2-1/apsi1)/phi^2
			varG2add3<-2*m^2*delta^2
			varG2<-varG2add1+varG2add2+varG2add3 # gamma_zero_teorico
			covG2null<- m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(1-exp(-delta*apsi1))*(exp(delta*apsi1)-1)/apsi1
			covG2<-function(h) covG2null*exp(-h*apsi1)
			CC<-matrix(numeric(q.back*q.back),q.back)
			
			CC<-covG2null*exp(-abs(row(CC)-col(CC))*delta*apsi1)
			
			diag(CC)<-numeric(q.back)+varG2
			b<-covG2null*exp(-(1:q.back)*delta*apsi1)
			a<-drop(solve(CC)%*%b)
			a0<-EG2*(1-sum(a))
			calcola.predittore <- function(x,sam) a0+a%*%sam[(x-1):(x-q.back)]
			#confronto.predittore<- function(x,sam) (sam[x]-(a0+a%*%sam[(x-1):(x-q.back)]))^2
			predittori<-sapply((q.back+1):length(samp2),calcola.predittore,sam=samp2)
			#c(sum((samp2[(q.back+1):length(samp2)]-predittori)^2),predittori)
			sum((samp2[(q.back+1):length(samp2)]-predittori)^2)
		}
	else stop("q.back argument needs to be smaller than the samp size")
	}
	else NA
}
)



setGeneric("MSPEestimation",function(init,samp,delta,q.back,levy) {
	init2<-standardGeneric("MSPEestimation")
	optimal<-optim(init2, MSPEcontrast, samp=samp, delta=delta, q.back=q.back, levy=levy)
	list(theta=COGprm(eta=optimal$par[2], beta=optimal$par[1], phi=optimal$par[3]),aux=optimal)
	}
)
setMethod("MSPEestimation",signature(init="vector"), function(init,samp,delta,q.back,levy) init)
setMethod("MSPEestimation",signature(init="COGprm"), function(init,samp,delta,q.back,levy) unname(as.vector(init)))
