#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// conventions
// I use the suffix "_param" for parameters passed from R
// They always are of type int* or double*
// Internal floating point variables and parameters are of type scalar, which should probably be double, but could be changed to long double if we want to


#include "coefficientFunctions.h"
#include "common.c"

void computeM0cases(
		const int *qParam,            // q
		const double *hParam,         // 
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		const double b[],              // 
		const double *azero_param,     //
		double M_result[]              // the resulting (q + 1) square matrix
		) {
	// assign parameters to variables
	const int q = *qParam;

	scalar *MM = calloc((q + 1)*(q + 1), sizeof(scalar));
	
        #define MM(j, k) MM[j*(q + 1) + k]
	const scalar azero = *azero_param;
        #include "coefficientsAsVariables.part.c"
	#include "casesM0.c"

	for (int jk = 0; jk < (q + 1)*(q + 1); jk++) {
		M_result[jk] = (double)MM[jk];
	}

	free(MM);
}
/*
void computeMINFcases(
		const int *qParam,            // q
		const double *hParam,         // 
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		const double b[],              // 
		const double *azero_param,     //
		double M_result[]              // the resulting (q + 1) square matrix
		) {
	// assign parameters to variables
	const int q = *qParam;

	scalar *MM = calloc((q + 1)*(q + 1), sizeof(scalar));
	
        #define MM(j, k) MM[j*(q + 1) + k]
	const scalar azero = *azero_param;
        #include "coefficientsAsVariables.part.c"
	#include "casesMINF.c"

	for (int jk = 0; jk < (q + 1)*(q + 1); jk++) {
		M_result[jk] = (double)MM[jk];
	}

	free(MM);
}
*/
void computeMNcases(
		const int *qParam,            // q
		const int *nParam,            // n
		const double *hParam,         // 
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		const double b[],              // 
		const double *azero_param,     //
		double M_result[]              // the resulting (q + 1) square matrix
		) {
	// assign parameters to variables
	const int q = *qParam;
	const int n = *nParam;
	if (n < q + 1) {
		printf("computeMNcases: parameter n = %d is smaller than (q + 1 = %d)\n", n, q + 1);
		return;
	}

	scalar *MM = calloc((q + 1)*(q + 1), sizeof(scalar));
	
        #define MM(j, k) MM[j*(q + 1) + k]
	const scalar azero = *azero_param;
        #include "coefficientsAsVariables.part.c"
	#include "casesMN.c"

	for (int jk = 0; jk < (q + 1)*(q + 1); jk++) {
		M_result[jk] = (double)MM[jk];
	}

	free(MM);
}

