#!/usr/bin/R

moment4 <- function(
		delta,r,t,
		beta, eta, tau, phi, 
		Psi1, Psi2, Psi3, Psi4) { 
	.C("moment4_R",
			as.double(delta),
			as.double(r), 
			as.double(t), 
			as.double(c(beta, eta, tau, phi)), 
			as.double(c(Psi1, Psi2, Psi3, Psi4)), 
			m = double(1)
	  )$m
}

moment6 <- function(
		delta, r, t,
		beta, eta, tau, phi, 
		Psi1, Psi2, Psi3, Psi4) {
	.C("moment6_R",
			as.double(delta),
			as.double(r), 
			as.double(t), 
			as.double(c(beta, eta, tau, phi)), 
			as.double(c(Psi1, Psi2, Psi3, Psi4)), 
			m = double(1)
	  )$m
}

moment8 <- function(
		delta,r,t,
		beta, eta, tau, phi, 
		Psi1, Psi2, Psi3, Psi4 
		) {
	.C("moment8_R",
			as.double(delta),
			as.double(r), 
			as.double(t), 
			as.double(c(beta, eta, tau, phi)), 
			as.double(c(Psi1, Psi2, Psi3, Psi4)), 
			m = double(1)
	  )$m
}

computeMklu <- function(d, n, h, 
        beta, eta, tau, phi, 
        Psi1, Psi2, Psi3, Psi4) {
    array(.C(
                "computeMklu", 
                as.integer(d), 
                as.integer(n), 
                as.double(h), 
                as.double(c(beta, eta, tau, phi)), 
                as.double(c(Psi1, Psi2, Psi3, Psi4)), 
                M=double((d + 2)^2)
            )$M, c(d + 2, d + 2))
}

computeM0 <- function(h, 
        beta, eta, tau, phi, 
        Psi1, Psi2, Psi3, Psi4, 
        a, a0) {
    q <- length(a);
    array(.C(
                "computeM0cases", 
                as.integer(q), 
                as.double(h), 
                as.double(c(beta, eta, tau, phi)), 
                as.double(c(Psi1, Psi2, Psi3, Psi4)), 
                as.double(c(-1, a)), 
                as.double(a0), 
                M=double((q+1)^2)
            )$M, c(q + 1, q + 1))
}

computeMN <- function(n, h, 
        beta, eta, tau, phi, 
        Psi1, Psi2, Psi3, Psi4, 
        a, a0) {
    q <- length(a);
    if (n <= 2*q + 2 ) {
        print("The n must be greater than 2q+1.");
        return;
    }
    m0 <- computeM0(h, beta, eta, tau, phi, Psi1, Psi2, Psi3, Psi4, a, a0);
    mn.asym <- array(.C(
                "computeMNcases", 
                as.integer(q), 
                as.integer(n), 
                as.double(h), 
                as.double(c(beta, eta, tau, phi)), 
                as.double(c(Psi1, Psi2, Psi3, Psi4)), 
                as.double(c(-1,a)), 
                as.double(a0), 
                M=double((q+1)^2)
                )$M, c(q + 1, q + 1));
    m0 + mn.asym + t(mn.asym)
}


computeMi <- function(q, n, i, h,
		beta, eta, tau, phi, 
		Psi1, Psi2, Psi3, Psi4, 
		b, azero) {
# put some consistency checks?
	Mi <- .C(
			"computeMi_R", 
			as.integer(i), 
			as.integer(q), 
			as.integer(n), 
			as.double(h), 
			as.double(c(beta, eta, tau, phi)), 
			as.double(c(Psi1, Psi2, Psi3, Psi4)), 
			as.double(b), 
			as.double(azero), 
			Mi=double((q+1)^2)
		  )$Mi;
	array(Mi, c(q + 1, q + 1))
}


computeMloop <- function(iter = 0, 
		q, n, h, 
		beta, eta, tau, phi, 
		Psi1, Psi2, Psi3, Psi4, 
		b, a0) {
	cat("Computing M[0]... ")
	Macc <- computeMi(q, n, 0, h, beta, eta, tau, phi, 
			Psi1, Psi2, Psi3, Psi4, b, a0)
	cat("done.\n")
	if (iter < 0) {
		cat("Argument iter is negative, setting it to 0\n")
		iter <- 0
	}
	else if (iter > n - q - 1) {
		cat("Argument iter is too big, setting it to", n - q - 1,"\n")
		iter <- n - q - 1
	}
	if (iter > 0) {
		for (i in 1:iter) {
			cat(sep = "", "Computing M[", i, "]... ")
			Mi <- computeMi(q, n, i, h, 
					beta, eta, tau, phi, 
					Psi1, Psi2, Psi3, Psi4, 
					b, a0)
			cat("done.\n")
			Macc <- Macc + Mi + t(Mi)
		}
	}
	Macc
}

