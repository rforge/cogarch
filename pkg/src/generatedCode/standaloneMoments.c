#include "common.c"
#include "coefficientFunctions.h"

void EG2(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	*result = ComputeEG2c(hParam,thetaParam,psiParam);
}

void EG2G2(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar delta1 = *delta1Param;
	
	*result = ComputeEG2G2c(hParam,thetaParam,psiParam) + ComputeEG2G2cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1);
}

void EG4(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	*result = ComputeEG4c(hParam,thetaParam,psiParam);
}

void EG2G2G2(const scalar *delta1Param, const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar delta1 = *delta1Param;
	const scalar delta2 = *delta2Param;
	
	*result = ComputeEG2G2G2c(hParam,thetaParam,psiParam) + ComputeEG2G2G2cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1) + ComputeEG2G2G2cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG2G2G2cAWAD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi1) + ComputeEG2G2G2cAWBD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi2);
}

void EG4G2(const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar delta2 = *delta2Param;
	
	*result = ComputeEG4G2c(hParam,thetaParam,psiParam) + ComputeEG4G2cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG4G2cBD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi2);
}

void EG2G4(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar delta1 = *delta1Param;
	
	*result = ComputeEG2G4c(hParam,thetaParam,psiParam) + ComputeEG2G4cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1);
}

void EG6(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	*result = ComputeEG6c(hParam,thetaParam,psiParam);
}

void EG2G2G2G2(const scalar *delta1Param, const scalar *delta2Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar Psi3 = psiParam[2];
	const scalar delta1 = *delta1Param;
	const scalar delta2 = *delta2Param;
	const scalar delta3 = *delta3Param;
	
	*result = ComputeEG2G2G2G2c(hParam,thetaParam,psiParam) + ComputeEG2G2G2G2cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1) + ComputeEG2G2G2G2cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG2G2G2G2cAH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi1) + ComputeEG2G2G2G2cAWAD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi1) + ComputeEG2G2G2G2cAWAH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta3*Psi1) + ComputeEG2G2G2G2cADAH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1 + delta3*Psi1) + ComputeEG2G2G2G2cAWADAH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi1 + delta3*Psi1) + ComputeEG2G2G2G2cAWBD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi2) + ComputeEG2G2G2G2cAWAHBD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta3*Psi1 + delta2*Psi2) + ComputeEG2G2G2G2cADBH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1 + delta3*Psi2) + ComputeEG2G2G2G2cAWADBH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi1 + delta3*Psi2) + ComputeEG2G2G2G2cAWBDBH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi2 + delta3*Psi2) + ComputeEG2G2G2G2cAWBDCH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi2 + delta3*Psi3);
}

void EG4G2G2(const scalar *delta2Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar Psi3 = psiParam[2];
	const scalar delta2 = *delta2Param;
	const scalar delta3 = *delta3Param;
	
	*result = ComputeEG4G2G2c(hParam,thetaParam,psiParam) + ComputeEG4G2G2cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG4G2G2cAH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi1) + ComputeEG4G2G2cADAH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1 + delta3*Psi1) + ComputeEG4G2G2cBD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi2) + ComputeEG4G2G2cAHBD(hParam,thetaParam,psiParam)*Exponential(delta3*Psi1 + delta2*Psi2) + ComputeEG4G2G2cADBH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1 + delta3*Psi2) + ComputeEG4G2G2cBDBH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi2 + delta3*Psi2) + ComputeEG4G2G2cBDCH(hParam,thetaParam,psiParam)*Exponential(delta2*Psi2 + delta3*Psi3);
}

void EG2G4G2(const scalar *delta1Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar Psi3 = psiParam[2];
	const scalar delta1 = *delta1Param;
	const scalar delta3 = *delta3Param;
	
	*result = ComputeEG2G4G2c(hParam,thetaParam,psiParam) + ComputeEG2G4G2cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1) + ComputeEG2G4G2cAH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi1) + ComputeEG2G4G2cAWAH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta3*Psi1) + ComputeEG2G4G2cBH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi2) + ComputeEG2G4G2cAWBH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta3*Psi2) + ComputeEG2G4G2cAWCH(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta3*Psi3);
}

void EG6G2(const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar Psi3 = psiParam[2];
	const scalar delta3 = *delta3Param;
	
	*result = ComputeEG6G2c(hParam,thetaParam,psiParam) + ComputeEG6G2cAH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi1) + ComputeEG6G2cBH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi2) + ComputeEG6G2cCH(hParam,thetaParam,psiParam)*Exponential(delta3*Psi3);
}

void EG2G2G4(const scalar *delta1Param, const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar delta1 = *delta1Param;
	const scalar delta2 = *delta2Param;
	
	*result = ComputeEG2G2G4c(hParam,thetaParam,psiParam) + ComputeEG2G2G4cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1) + ComputeEG2G2G4cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG2G2G4cAWAD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi1) + ComputeEG2G2G4cAWBD(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1 + delta2*Psi2);
}

void EG4G4(const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar Psi2 = psiParam[1];
	const scalar delta2 = *delta2Param;
	
	*result = ComputeEG4G4c(hParam,thetaParam,psiParam) + ComputeEG4G4cAD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi1) + ComputeEG4G4cBD(hParam,thetaParam,psiParam)*Exponential(delta2*Psi2);
}

void EG2G6(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	const scalar Psi1 = psiParam[0];
	const scalar delta1 = *delta1Param;
	
	*result = ComputeEG2G6c(hParam,thetaParam,psiParam) + ComputeEG2G6cAW(hParam,thetaParam,psiParam)*Exponential(delta1*Psi1);
}

void EG8(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result)
{
	*result = ComputeEG8c(hParam,thetaParam,psiParam);
}
