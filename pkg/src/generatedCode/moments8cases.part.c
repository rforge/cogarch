case 0: // EG2G2G2G2
return coeff[EG2G2G2G2c] + coeff[EG2G2G2G2cAW]*psiexp[3*delta1] + coeff[EG2G2G2G2cAD]*psiexp[3*delta2] + coeff[EG2G2G2G2cAWAD]*psiexp[3*delta1]*psiexp[3*delta2] + coeff[EG2G2G2G2cAWBD]*psiexp[3*delta1]*psiexp[1 + 3*delta2] + coeff[EG2G2G2G2cAH]*psiexp[3*delta3] + coeff[EG2G2G2G2cAWAH]*psiexp[3*delta1]*psiexp[3*delta3] + coeff[EG2G2G2G2cADAH]*psiexp[3*delta2]*psiexp[3*delta3] + coeff[EG2G2G2G2cAWADAH]*psiexp[3*delta1]*psiexp[3*delta2]*psiexp[3*delta3] + coeff[EG2G2G2G2cAWAHBD]*psiexp[3*delta1]*psiexp[1 + 3*delta2]*psiexp[3*delta3] + coeff[EG2G2G2G2cADBH]*psiexp[3*delta2]*psiexp[1 + 3*delta3] + coeff[EG2G2G2G2cAWADBH]*psiexp[3*delta1]*psiexp[3*delta2]*psiexp[1 + 3*delta3] + coeff[EG2G2G2G2cAWBDBH]*psiexp[3*delta1]*psiexp[1 + 3*delta2]*psiexp[1 + 3*delta3] + coeff[EG2G2G2G2cAWBDCH]*psiexp[3*delta1]*psiexp[1 + 3*delta2]*psiexp[2 + 3*delta3];
case 1: // EG4G2G2
return coeff[EG4G2G2c] + coeff[EG4G2G2cAD]*psiexp[3*delta2] + coeff[EG4G2G2cBD]*psiexp[1 + 3*delta2] + coeff[EG4G2G2cAH]*psiexp[3*delta3] + coeff[EG4G2G2cADAH]*psiexp[3*delta2]*psiexp[3*delta3] + coeff[EG4G2G2cAHBD]*psiexp[1 + 3*delta2]*psiexp[3*delta3] + coeff[EG4G2G2cADBH]*psiexp[3*delta2]*psiexp[1 + 3*delta3] + coeff[EG4G2G2cBDBH]*psiexp[1 + 3*delta2]*psiexp[1 + 3*delta3] + coeff[EG4G2G2cBDCH]*psiexp[1 + 3*delta2]*psiexp[2 + 3*delta3];
case 2: // EG2G4G2
return coeff[EG2G4G2c] + coeff[EG2G4G2cAW]*psiexp[3*delta1] + coeff[EG2G4G2cAH]*psiexp[3*delta3] + coeff[EG2G4G2cAWAH]*psiexp[3*delta1]*psiexp[3*delta3] + coeff[EG2G4G2cBH]*psiexp[1 + 3*delta3] + coeff[EG2G4G2cAWBH]*psiexp[3*delta1]*psiexp[1 + 3*delta3] + coeff[EG2G4G2cAWCH]*psiexp[3*delta1]*psiexp[2 + 3*delta3];
case 3: // EG6G2
return coeff[EG6G2c] + coeff[EG6G2cAH]*psiexp[3*delta3] + coeff[EG6G2cBH]*psiexp[1 + 3*delta3] + coeff[EG6G2cCH]*psiexp[2 + 3*delta3];
case 4: // EG2G2G4
return coeff[EG2G2G4c] + coeff[EG2G2G4cAW]*psiexp[3*delta1] + coeff[EG2G2G4cAD]*psiexp[3*delta2] + coeff[EG2G2G4cAWAD]*psiexp[3*delta1]*psiexp[3*delta2] + coeff[EG2G2G4cAWBD]*psiexp[3*delta1]*psiexp[1 + 3*delta2];
case 5: // EG4G4
return coeff[EG4G4c] + coeff[EG4G4cAD]*psiexp[3*delta2] + coeff[EG4G4cBD]*psiexp[1 + 3*delta2];
case 6: // EG2G6
return coeff[EG2G6c] + coeff[EG2G6cAW]*psiexp[3*delta1];
case 7: // EG8
return coeff[EG8c];
