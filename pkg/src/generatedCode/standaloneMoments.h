#include "common.h"

void EG2(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G2(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG4(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G2G2(const scalar *delta1Param, const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG4G2(const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G4(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG6(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G2G2G2(const scalar *delta1Param, const scalar *delta2Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG4G2G2(const scalar *delta2Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G4G2(const scalar *delta1Param, const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG6G2(const scalar *delta3Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G2G4(const scalar *delta1Param, const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG4G4(const scalar *delta2Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG2G6(const scalar *delta1Param, const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
void EG8(const scalar *hParam, const scalar *thetaParam, const scalar *psiParam, scalar *result);
