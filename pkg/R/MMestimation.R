
MMestimation<-function(x,d=sqrt(length(x)),explicit=TRUE){ 
	mu<-mean(x^2)
	gamma0<-mean(x^4)-mu^2
	logs<-log(abs(acf(x^2,plot=FALSE,lag.max=d,type="correlation")$acf[-1])) #NB qui fare abs non mi piace!! per ora è così
	pesi<- (1:d)-(d+1)/2
	ml<-mean(logs)
	if(explicit){
		p1<- -sum((logs-ml)*pesi)/sum(pesi^2)
		k1<- exp(ml+((d+1)*p1)/2)
	}
	else{
		Hfun<-function(theta) sum((logs- log(theta[2])+theta[1]*(1:length(logs)))^2)
		sol<-optim(c(0.1,0.1),Hfun,method="L-BFGS-B",lower=c(0.0001,0.0001))
		p1<-sol$par[1]
		k1<-sol$par[2]
	}
	posp<- p1>0
	M1<-gamma0-2*mu^2-6*k1*gamma0*(1-p1-exp(-p1))/((1-exp(p1))*(1-exp(-p1)))
	M2<-2*k1*gamma0*p1/(M1*(exp(p1)-1)*(1-exp(-p1)))
	beta<-p1*mu
	eta<-p1*sqrt(1+M2)
	phi<-eta-p1
	theta<-COGprm(eta=eta, beta=beta, phi=phi)
	return(list(theta=theta,aux=posp))
}

