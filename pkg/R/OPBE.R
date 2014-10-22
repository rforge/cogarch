setGeneric("OPBEF", function(theta, samp, delta, q.back, ...){
standardGeneric("OPBEF")
}
)

setMethod("OPBEF",c("COGprm", "vector", "numeric", "numeric"), function(theta, samp, delta, q.back,...){#... dovrebbe per ora contenere il vg o altra specifica del levy 
#delta is related to the fact that observations are from G^(delta)_t=G_t-G_(t-delta) with G_t solving dG_t=sigma_t dL_t #I will in part mimic notations in Sorensen's paper
	if (theta@beta>0.000001 && theta@eta>0.000001 &&theta@phi>0.000001){
		if(q.back<length(samp)){
			beta<-theta@beta
			eta<-theta@eta
			phi<-theta@phi
			lista<-setPsi(theta,...)
			psi1<-lista$psi1
			psi2<-lista$psi2
			psi3<-lista$psi3
			psi4<-lista$psi4
			tau<-lista$tau
			samp2<-samp^2 #in input vogliamo i log-returns, non i quadrati, ma se si hanno i log-prezzi bisogna già applicare le diff prima
			lung<-length(samp)
			if (psi4>=0) return(NA)
			apsi1<-abs(psi1)
			apsi2<-abs(psi2)
			EG2<-beta*delta/apsi1
			m<-beta/apsi1
			v<-2*beta^2/abs(psi1*psi2)
			varG2add1<- 6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)
			varG2add2<- 2*delta*beta^2*(2/apsi2-1/apsi1)/phi^2
			varG2add3<-2*m^2*delta^2
			varG2<-varG2add1+varG2add2+varG2add3 
			covG2null<- m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(1-exp(-delta*apsi1))*(exp(delta*apsi1)-1)/apsi1
			covG2<-function(h) covG2null*exp(-h*apsi1)
			CC<-matrix(numeric(q.back*q.back),q.back)
			CCmat_pos<- -abs(row(CC)-col(CC))
			CCmat_exp<-exp(CCmat_pos*delta*apsi1)
			CC<-covG2null*CCmat_exp
			diag(CC)<-numeric(q.back)+varG2
			b<-covG2null*exp(-(1:q.back)*delta*apsi1)
			CCinv<-try(solve(CC))
			if(!all(is.finite(CCinv))) {print(paste(theta,"CC")); return(NA)}
			a<-CCinv%*%b
			a0<-EG2*(1-sum(a))
			tildeC<-rbind(EG2,cbind(EG2,CC+EG2^2)) #Ctilte prima 3.13 S.
			tildeC[1,1]<-1
			
			#le derivate sono cacolate assumendo gli psi negativi cioè usando apsi=-psi
			apsi1_eta<- 1
			apsi2_eta<- 2
			apsi1_phi<- (apsi1-eta)/phi
			apsi2_phi<- 2*(apsi2-2*apsi1)/phi + 2*apsi1_phi
			
			EG2_beta<- EG2/beta
			EG2_eta<- (-EG2/apsi1)*apsi1_eta
			EG2_phi<- (-EG2/apsi1)*apsi1_phi
			covG2null_beta<- (2/beta)*covG2null
			covG2null_apsi1<- covG2null*(-3/apsi1+ 1/((2/apsi2-1/apsi1)*apsi1^2)+ 2*delta*exp(-delta*apsi1)/(1-exp(-delta*apsi1))+delta+ 2/(2*apsi1+phi))
			covG2null_apsi2<- -2*covG2null/((2/apsi2-1/apsi1)*apsi2^2) ###NB a me sembra corretto

			covG2null_eta<- covG2null_apsi1*apsi1_eta+covG2null_apsi2*apsi2_eta
			covG2null_phi<- -2*apsi1*covG2null/(phi*(2*apsi1+phi))+ covG2null_apsi1*apsi1_phi +covG2null_apsi2*apsi2_phi
			varG2_beta<-2*varG2/beta
			varG2_apsi1<-varG2add1*(-3/apsi1+ 1/((2/apsi2-1/apsi1)*apsi1^2)+ delta*(1-exp(-delta*apsi1))/(delta*apsi1-1+exp(-delta*apsi1))+ 2/(2*apsi1+phi)) +2*delta*beta^2/(apsi1^2*phi^2) - 4*beta^2*delta^2/apsi1^3
			varG2_apsi2<- -2*(varG2add1+varG2add2)/((2/apsi2-1/apsi1)*apsi2^2)
			varG2_eta<- varG2_apsi1*apsi1_eta+ varG2_apsi2*apsi2_eta
			varG2_phi<- -2*apsi1*varG2add1/(phi*(2*apsi1+phi)) -2* varG2add2/phi + varG2_apsi1*apsi1_phi + varG2_apsi2*apsi2_phi			
			
			CCmat_exp_apsi<-exp(CCmat_pos*delta*apsi1)*CCmat_pos*delta
			CCmat_exp_eta<-CCmat_exp_apsi*apsi1_eta
			CCmat_exp_phi<-CCmat_exp_apsi*apsi1_phi
			
			CC_beta<-covG2null_beta*CCmat_exp
			diag(CC_beta)<-numeric(q.back)+varG2_beta
			CC_eta<-covG2null_eta*CCmat_exp+ covG2null*CCmat_exp_eta
			diag(CC_eta)<-numeric(q.back)+varG2_eta
			CC_phi<-covG2null_phi*CCmat_exp+ covG2null*CCmat_exp_phi
			diag(CC_phi)<-numeric(q.back)+varG2_phi
			b_beta<- covG2null_beta*exp(-(1:q.back)*delta*apsi1)
			b_eta<- covG2null_eta*exp(-(1:q.back)*delta*apsi1)+covG2null*exp(-(1:q.back)*delta*apsi1)*(-(1:q.back)*delta)*apsi1_eta
			b_phi<- covG2null_phi*exp(-(1:q.back)*delta*apsi1)+covG2null*exp(-(1:q.back)*delta*apsi1)*(-(1:q.back)*delta)*apsi1_phi
			a_grad_zero_escluso <- CCinv%*%cbind(b_beta- CC_beta%*%a, b_eta- CC_eta%*%a, b_phi- CC_phi%*%a) #matrice qxp
			a_beta <- CCinv%*%(b_beta- CC_beta%*%a)
			a_eta <- CCinv%*%(b_eta- CC_eta%*%a)
			a_phi <- CCinv%*%(b_phi- CC_phi%*%a)
			a0_beta<-EG2_beta*(1-sum(a))-EG2*sum(a_beta)
			a0_eta<-EG2_eta*(1-sum(a))-EG2*sum(a_eta)
			a0_phi<-EG2_phi*(1-sum(a))-EG2*sum(a_phi)
			a_grad <- rbind(c(a0_beta,a0_eta, a0_phi),a_grad_zero_escluso)# matrice (q+1)xp
			sseq<- -1:-q.back
			traspa<-t(a)
			somma.H<-0
			somma0<-0
			for(i in (q.back+1):length(samp2)){
				predittore<-a0+ traspa%*%samp2[i+sseq]
				somma0<-somma0+(samp2[i]-predittore)
				somma.H<-somma.H+samp2[i+sseq]*(samp2[i]-predittore)
			}
			somma.H<-c(somma0,somma.H)
			#U dell 3.12 S ed è una matrice (q+1)xp
			U <- tildeC %*% a_grad
			M0<-computeMi(q.back, lung, 0, 1, beta, eta, tau, phi, psi1, psi2, psi3, psi4, c(-1,a), a0)
			M0inv<-try(solve(M0))
			if(!all(is.finite(M0inv))) {print(theta); return(NA)}
			A<- t(U)%*%M0inv
			return(A%*%somma.H)
			#estimating<- A%*%somma.H
			#return(t(estimating)%*%estimating)
		}
		else stop("q.back argument needs to be smaller than the sample size")
	} else {return(NA)}  
}
)

setMethod("OPBEF",c("vector", "vector", "numeric", "numeric"), function(theta, samp, delta, q.back,...){#... dovrebbe per ora contenere il vg o altra specifica del levy 
#delta is related to the fact that observations are from G^(delta)_t=G_t-G_(t-delta) with G_t solving dG_t=sigma_t dL_t #I will in part mimic notations in Sorensen's paper
	if (theta[1]>0.000001 && theta[2]>0.000001 &&theta[3]>0.000001){
		if(q.back<length(samp)){
			beta<-theta[1]
			eta<-theta[2]
			phi<-theta[3]
			lista<-setPsi(COGprm(eta=eta, beta=beta, phi=phi),...)
			psi1<-lista$psi1
			psi2<-lista$psi2
			psi3<-lista$psi3
			psi4<-lista$psi4
			tau<-lista$tau
			samp2<-samp^2 #in input vogliamo i log-returns, non i quadrati, ma se si hanno i log-prezzi bisogna già applicare le diff prima
			lung<-length(samp)
			if (psi4>=0) return(NA)
			apsi1<-abs(psi1)
			apsi2<-abs(psi2)
			EG2<-beta*delta/apsi1
			m<-beta/apsi1
			v<-2*beta^2/abs(psi1*psi2)
			varG2add1<- 6*m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(delta-(1-exp(-delta*apsi1))/apsi1)
			varG2add2<- 2*delta*beta^2*(2/apsi2-1/apsi1)/phi^2
			varG2add3<-2*m^2*delta^2
			varG2<-varG2add1+varG2add2+varG2add3 
			covG2null<- m^2*((2*apsi1+phi)/phi)*(2/apsi2-1/apsi1)*(1-exp(-delta*apsi1))*(exp(delta*apsi1)-1)/apsi1
			covG2<-function(h) covG2null*exp(-h*apsi1)
			CC<-matrix(numeric(q.back*q.back),q.back)
			CCmat_pos<- -abs(row(CC)-col(CC))
			CCmat_exp<-exp(CCmat_pos*delta*apsi1)
			CC<-covG2null*CCmat_exp
			diag(CC)<-numeric(q.back)+varG2
			b<-covG2null*exp(-(1:q.back)*delta*apsi1)
			CCinv<-try(solve(CC))
			if(!all(is.finite(CCinv))) {print(paste(theta,"CC")); return(NA)}
			a<-CCinv%*%b
			a0<-EG2*(1-sum(a))
			tildeC<-rbind(EG2,cbind(EG2,CC+EG2^2)) #Ctilte prima 3.13 S.
			tildeC[1,1]<-1
			
			#le derivate sono cacolate assumendo gli psi negativi cioè usando apsi=-psi
			apsi1_eta<- 1
			apsi2_eta<- 2
			apsi1_phi<- (apsi1-eta)/phi
			apsi2_phi<- 2*(apsi2-2*apsi1)/phi + 2*apsi1_phi
			
			EG2_beta<- EG2/beta
			EG2_eta<- (-EG2/apsi1)*apsi1_eta
			EG2_phi<- (-EG2/apsi1)*apsi1_phi
			covG2null_beta<- (2/beta)*covG2null
			covG2null_apsi1<- covG2null*(-3/apsi1+ 1/((2/apsi2-1/apsi1)*apsi1^2)+ 2*delta*exp(-delta*apsi1)/(1-exp(-delta*apsi1))+delta+ 2/(2*apsi1+phi))
			covG2null_apsi2<- -2*covG2null/((2/apsi2-1/apsi1)*apsi2^2) ###NB a me sembra corretto

			covG2null_eta<- covG2null_apsi1*apsi1_eta+covG2null_apsi2*apsi2_eta
			covG2null_phi<- -2*apsi1*covG2null/(phi*(2*apsi1+phi))+ covG2null_apsi1*apsi1_phi +covG2null_apsi2*apsi2_phi
			varG2_beta<-2*varG2/beta
			varG2_apsi1<-varG2add1*(-3/apsi1+ 1/((2/apsi2-1/apsi1)*apsi1^2)+ delta*(1-exp(-delta*apsi1))/(delta*apsi1-1+exp(-delta*apsi1))+ 2/(2*apsi1+phi)) +2*delta*beta^2/(apsi1^2*phi^2) - 4*beta^2*delta^2/apsi1^3
			varG2_apsi2<- -2*(varG2add1+varG2add2)/((2/apsi2-1/apsi1)*apsi2^2)
			varG2_eta<- varG2_apsi1*apsi1_eta+ varG2_apsi2*apsi2_eta
			varG2_phi<- -2*apsi1*varG2add1/(phi*(2*apsi1+phi)) -2* varG2add2/phi + varG2_apsi1*apsi1_phi + varG2_apsi2*apsi2_phi			
			
			CCmat_exp_apsi<-exp(CCmat_pos*delta*apsi1)*CCmat_pos*delta
			CCmat_exp_eta<-CCmat_exp_apsi*apsi1_eta
			CCmat_exp_phi<-CCmat_exp_apsi*apsi1_phi
			
			CC_beta<-covG2null_beta*CCmat_exp
			diag(CC_beta)<-numeric(q.back)+varG2_beta
			CC_eta<-covG2null_eta*CCmat_exp+ covG2null*CCmat_exp_eta
			diag(CC_eta)<-numeric(q.back)+varG2_eta
			CC_phi<-covG2null_phi*CCmat_exp+ covG2null*CCmat_exp_phi
			diag(CC_phi)<-numeric(q.back)+varG2_phi
			b_beta<- covG2null_beta*exp(-(1:q.back)*delta*apsi1)
			b_eta<- covG2null_eta*exp(-(1:q.back)*delta*apsi1)+covG2null*exp(-(1:q.back)*delta*apsi1)*(-(1:q.back)*delta)*apsi1_eta
			b_phi<- covG2null_phi*exp(-(1:q.back)*delta*apsi1)+covG2null*exp(-(1:q.back)*delta*apsi1)*(-(1:q.back)*delta)*apsi1_phi
			a_grad_zero_escluso <- CCinv%*%cbind(b_beta- CC_beta%*%a, b_eta- CC_eta%*%a, b_phi- CC_phi%*%a) #matrice qxp
			a_beta <- CCinv%*%(b_beta- CC_beta%*%a)
			a_eta <- CCinv%*%(b_eta- CC_eta%*%a)
			a_phi <- CCinv%*%(b_phi- CC_phi%*%a)
			a0_beta<-EG2_beta*(1-sum(a))-EG2*sum(a_beta)
			a0_eta<-EG2_eta*(1-sum(a))-EG2*sum(a_eta)
			a0_phi<-EG2_phi*(1-sum(a))-EG2*sum(a_phi)
			a_grad <- rbind(c(a0_beta,a0_eta, a0_phi),a_grad_zero_escluso)# matrice (q+1)xp
			sseq<- -1:-q.back
			traspa<-t(a)
			somma.H<-0
			somma0<-0
			for(i in (q.back+1):length(samp2)){
				predittore<-a0+ traspa%*%samp2[i+sseq]
				somma0<-somma0+(samp2[i]-predittore)
				somma.H<-somma.H+samp2[i+sseq]*(samp2[i]-predittore)
			}
			somma.H<-c(somma0,somma.H)
			#U dell 3.12 S ed è una matrice (q+1)xp
			U <- tildeC %*% a_grad
			M0<-computeMi(q.back, lung, 0, 1, beta, eta, tau, phi, psi1, psi2, psi3, psi4, c(-1,a), a0)
			M0inv<-try(solve(M0))
			if(!all(is.finite(M0inv))) {print(theta); return(NA)}
			A<- t(U)%*%M0inv
			return(A%*%somma.H)
			#estimating<- A%*%somma.H
			#return(t(estimating)%*%estimating)
		}
		else stop("q.back argument needs to be smaller than the sample size")
	} else {return(NA)}  
}
)



OPBcontrast<-function(theta, samp, delta, q.back, levy) {
	estimating<-OPBEF(theta, samp, delta, q.back, levy)
	t(estimating)%*%estimating
	}

	
setGeneric("OPBestimation",function(init,samp,delta,q.back,levy) {
	init2<-standardGeneric("OPBestimation")
	optimal<-optim(init2, OPBcontrast, samp=samp, delta=delta, q.back=q.back, levy=levy)
	list(theta=COGprm(eta=optimal$par[2], beta=optimal$par[1], phi=optimal$par[3]),aux=optimal)
	}
)
setMethod("OPBestimation",signature(init="vector"), function(init,samp,delta,q.back,levy) init)
setMethod("OPBestimation",signature(init="COGprm"), function(init,samp,delta,q.back,levy) unname(as.vector(init)))



####################################################











