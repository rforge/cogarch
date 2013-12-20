case 0: // EG2G2G2
return coeff[EG2G2G2c] + coeff[EG2G2G2cAW]*psiexp[3*delta1] + coeff[EG2G2G2cAD]*psiexp[3*delta2] + coeff[EG2G2G2cAWAD]*psiexp[3*delta1]*psiexp[3*delta2] + coeff[EG2G2G2cAWBD]*psiexp[3*delta1]*psiexp[1 + 3*delta2];
case 1: // EG4G2
return coeff[EG4G2c] + coeff[EG4G2cAD]*psiexp[3*delta2] + coeff[EG4G2cBD]*psiexp[1 + 3*delta2];
case 2: // EG2G4
return coeff[EG2G4c] + coeff[EG2G4cAW]*psiexp[3*delta1];
case 3: // EG6
return coeff[EG6c];
