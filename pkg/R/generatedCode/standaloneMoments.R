EG2 <- function(h, theta, psi)
.C("EG2", as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G2 <- function(delta1, h, theta, psi)
.C("EG2G2", as.double(delta1), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG4 <- function(h, theta, psi)
.C("EG4", as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G2G2 <- function(delta1, delta2, h, theta, psi)
.C("EG2G2G2", as.double(delta1), as.double(delta2), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG4G2 <- function(delta2, h, theta, psi)
.C("EG4G2", as.double(delta2), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G4 <- function(delta1, h, theta, psi)
.C("EG2G4", as.double(delta1), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG6 <- function(h, theta, psi)
.C("EG6", as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G2G2G2 <- function(delta1, delta2, delta3, h, theta, psi)
.C("EG2G2G2G2", as.double(delta1), as.double(delta2), as.double(delta3), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG4G2G2 <- function(delta2, delta3, h, theta, psi)
.C("EG4G2G2", as.double(delta2), as.double(delta3), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G4G2 <- function(delta1, delta3, h, theta, psi)
.C("EG2G4G2", as.double(delta1), as.double(delta3), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG6G2 <- function(delta3, h, theta, psi)
.C("EG6G2", as.double(delta3), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G2G4 <- function(delta1, delta2, h, theta, psi)
.C("EG2G2G4", as.double(delta1), as.double(delta2), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG4G4 <- function(delta2, h, theta, psi)
.C("EG4G4", as.double(delta2), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG2G6 <- function(delta1, h, theta, psi)
.C("EG2G6", as.double(delta1), as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result

EG8 <- function(h, theta, psi)
.C("EG8", as.double(h), as.double(theta), as.double(psi), result=as.double(numeric(1)))$result
