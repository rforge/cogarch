setGeneric("PML", function(theta, samp, delta){
standardGeneric("PML")
}
)

setMethod("PML",c("COGprm", "vector", "numeric"), function(theta,samp,delta){   #NB x sono diff(G)  x2 sono x^2 
	#qui delta è costante. il metodo vale anche se non sono equispaziati, ma io non lo ho scritto in quest'ottica.
	######è corretta nel rho2, e dà output solo se eta<phi, per cui il processo è stazionario
	x<-samp
	beta<-theta@beta
	eta<-theta@eta
	phi<-theta@phi
	if(eta>phi) {
		x2<-x^2
		n<-length(x2)
		sigma2<-numeric(n) #NB l'indicizzazione non è corente con quella degli articoli: sigma[i] e x[i-1]si riferiscono allo stesso tempo reale
		#rho2<-numeric(n) #il tempo di rho invece è allineato a quello di x 
		sigma2[1]<-beta/(eta-phi) 
		#rho2[1]<- (sigma2[1] -beta/(eta-phi))*((exp(-(eta-phi)*delta)-1)/(-eta+phi))+beta*delta/(eta-phi)
		for (i in 2:n){
		sigma2[i]<- beta*delta +exp(-eta*delta)*sigma2[i-1]+phi*exp(-eta*delta)*x2[i-1]
		#rho2[i]<- (sigma2[i] -beta/(eta-phi))*((exp((eta-phi)*delta)-1)/(eta-phi))+beta*delta/(eta-phi)#NB può essere vettorizzato!
		}
		rho2<- (sigma2 -beta/(eta-phi))*((exp(-(eta-phi)*delta)-1)/(-eta+phi))+beta*delta/(eta-phi)
		return(sum(x2/rho2+log(rho2)))
	}

	else NA
}	
)	

setMethod("PML",c("vector", "vector", "numeric"), function(theta,samp,delta){   #NB x sono diff(G)  x2 sono x^2 
	#qui delta è costante. il metodo vale anche se non sono equispaziati, ma io non lo ho scritto in quest'ottica.
	######è corretta nel rho2, e dà output solo se eta<phi, per cui il processo è stazionario
	x<-samp
	beta<-theta[1]
	eta<-theta[2]
	phi<-theta[3]
	if(eta>phi) {
		x2<-x^2
		n<-length(x2)
		sigma2<-numeric(n) #NB l'indicizzazione non è corente con quella degli articoli: sigma[i] e x[i-1]si riferiscono allo stesso tempo reale
		#rho2<-numeric(n) #il tempo di rho invece è allineato a quello di x 
		sigma2[1]<-beta/(eta-phi) 
		#rho2[1]<- (sigma2[1] -beta/(eta-phi))*((exp(-(eta-phi)*delta)-1)/(-eta+phi))+beta*delta/(eta-phi)
		for (i in 2:n){
		sigma2[i]<- beta*delta +exp(-eta*delta)*sigma2[i-1]+phi*exp(-eta*delta)*x2[i-1]
		#rho2[i]<- (sigma2[i] -beta/(eta-phi))*((exp((eta-phi)*delta)-1)/(eta-phi))+beta*delta/(eta-phi)#NB può essere vettorizzato!
		}
		rho2<- (sigma2 -beta/(eta-phi))*((exp(-(eta-phi)*delta)-1)/(-eta+phi))+beta*delta/(eta-phi)
		return(sum(x2/rho2+log(rho2)))
	}

	else NA
}	
)	


# tutti i metodi di stima mangiano come campione i diff(G) in condizione stazionaria
#theta,samp,delta,q.back

setGeneric("PMLestimation",function(init,samp,delta) {
	init2<-standardGeneric("PMLestimation")
	optimal<-optim(init2, PML, samp=samp, delta=delta)
	list(theta=COGprm(eta=optimal$par[2], beta=optimal$par[1], phi=optimal$par[3]),aux=optimal)
	}
)
setMethod("PMLestimation",signature(init="vector"), function(init,samp,delta) init)
setMethod("PMLestimation",signature(init="COGprm"), function(init,samp,delta) unname(as.vector(init)))
