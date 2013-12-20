#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// conventions
// I use the suffix "Param" for parameters passed from R
// They always are of type int* or double*
// Internal floating point variables and parameters are of type scalar, which should probably be double, but could be changed to long double if we want to


//#include "coefficients.h"
#include "common.c"
#include "coefficientFunctions.h"
#include "coefficientsAsArray.h"

// Fast sorting for four integers.
// Since we only have to sort four elements we can use sorting networks which should be the best choice for our needs
// For more information see Knuth's Art of Computer Programming, Volume III

static inline void Sort2(int *p0, int *p1) {
#if SWAP == 0
#define min(x, y) (y ^ ((x ^ y) & -(x < y)))
#define max(x, y) (x ^ ((x ^ y) & -(x < y)))
	const int temp = min(*p0, *p1);
	*p1 = max(*p0, *p1);
	*p0 = temp;
#elif SWAP == 1
	if (*p1 < *p0) {
		int temp = *p0;
		*p1 = *p0;
		*p0 = temp;
	}
#else
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))
	const int temp = min(*p0, *p1);
	*p1 = max(*p0, *p1);
	*p0 = temp;

#endif
}

scalar moment4(const int delta, const scalar psiexp[], const scalar coeff[]) {
	const int delta1 = abs(delta);
	switch ((delta1 == 0) << 0) {
#include "moments4cases.part.c"
	}
	assert(0);
	return 0;
}

scalar moment6(int time[], const scalar psiexp[], const scalar coeff[]) {
	Sort2(time + 0, time + 1);
	Sort2(time + 1, time + 2);
	Sort2(time + 0, time + 1);
	const int delta1 = time[2] - time[1];
	const int delta2 = time[1] - time[0];
	assert(delta1 >= 0);
	assert(delta2 >= 0);
	switch (((delta1 == 0) << 0) | ((delta2 == 0) << 1)) {
#include "moments6cases.part.c"
	}
	assert(0);
	return 0;
}

scalar moment8(int time[], const scalar psiexp[], const scalar coeff[]) {
	Sort2(time + 0, time + 1);
	Sort2(time + 2, time + 3);
	Sort2(time + 0, time + 2);
	Sort2(time + 1, time + 3);
	Sort2(time + 1, time + 2);
	const int delta1 = time[3] - time[2];
	const int delta2 = time[2] - time[1];
	const int delta3 = time[1] - time[0];
	assert(delta1 >= 0);
	assert(delta2 >= 0);
	assert(delta3 >= 0);
	switch (((delta1 == 0) << 0) | ((delta2 == 0) << 1) | ((delta3 == 0) << 2)) {
#include "moments8cases.part.c"
	}
	assert(0);
	return 0;
}

#define moment4(delta) moment4(delta, psiexp, coeff)
#define moment6(time) moment6(time, psiexp, coeff)
#define moment8(time) moment8(time, psiexp, coeff)
/*
void moment4_R(
		const int *delta, 
		const double *rParam,         // r
		const double *tParam,         // t
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		double *moment_result
		) {
	// assign parameters to variables
	const int n = *delta + 1;

	// precompute exponentials and coefficients
	scalar *psiexp = computePsiExp(&n, psiParam);
	scalar *coeff = computeCoefficients(rParam, tParam, thetaParam, psiParam);

	*moment_result = moment4(*delta);

	free(psiexp);
	free(coeff);
}

void moment6_R(
		const int time[], 
		const int *nParam,
		const double *rParam,         // r
		const double *tParam,         // t
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		double *moment_result
		) {
	// assign parameters to variables
	int _time[3];
	_time[0] = time[0];
	_time[1] = time[1];
	_time[2] = time[2];

	// precompute exponentials and coefficients
	scalar *psiexp = computePsiExp(nParam, psiParam);
	scalar *coeff = computeCoefficients(rParam, tParam, thetaParam, psiParam);

	*moment_result = moment6(_time);

	free(psiexp);
	free(coeff);
}

void moment8_R(
		const int time[], 
		const int *nParam,
		const double *rParam,         // r
		const double *tParam,         // t
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		double *moment_result
		) {
	// precompute exponentials and coefficients
	scalar *psiexp = computePsiExp(nParam, psiParam);
	scalar *coeff = computeCoefficients(rParam, tParam, thetaParam, psiParam);

	int _time[4];
	_time[0] = time[0];
	_time[1] = time[1];
	_time[2] = time[2];
	_time[3] = time[3];

	*moment_result = moment8(_time);

	free(psiexp);
	free(coeff);
}

*/

// precomputation: exponentials of multiples of Psi1, Psi2, Psi3
// this is a very significant optimization, from 10 minutes to 34 seconds in the test case
scalar *computePsiExp(const int *nParam, const scalar psiParam[4]) {
	const int n = *nParam;
	const scalar Psi1 = psiParam[0]; 
	const scalar Psi2 = psiParam[1]; 
	const scalar Psi3 = psiParam[2];

	// the maximum delta seems to be (n - 1)
	scalar *psiexp = (scalar*)malloc(sizeof(scalar)*3*n);
	assert(psiexp != NULL);

	for (int delta = 0; delta < n; delta++) {
		psiexp[3*delta + 0] = Exponential(Psi1*delta); 
		psiexp[3*delta + 1] = Exponential(Psi2*delta); 
		psiexp[3*delta + 2] = Exponential(Psi3*delta); 
	}
	return psiexp;
}



void computeMi(const int i, const int q, const int n, const double b[], const scalar a_zero, const scalar psiexp[], const scalar coeff[], scalar Mi[]) {
#define Mi(j, k) Mi[j*(q + 1) + k]

	{// j = 0, k = 0

		scalar sum_kappa = 0;
		for (int kappa = 0; kappa <= q; kappa++) {
			sum_kappa += b[kappa];
		}

		scalar sum_nu_kappa = 0;
		for (int kappa = 0; kappa <= q; kappa++) {
			for (int nu = 0; nu <= q; nu++) {
				sum_nu_kappa += b[kappa]*b[nu]*moment4(i - kappa + nu);
			}
		}
		Mi(0, 0) = (n - q - i) * (sum_nu_kappa + a_zero*2*sum_kappa*coeff[EG2c] + a_zero*a_zero) / (n - q); // half of 3.14
	}

	// k = 0, j > 0
	for (int j = 1; j <= q; j++) {
		scalar sum_nu_kappa = 0;
		scalar sum_kappa = 0;
		for (int kappa = 0; kappa <= q; kappa++) {
			for (int nu = 0; nu <= q; nu++) {
				int times[3];
				times[0] = -j;
				times[1] = -nu;
				times[2] = i - kappa;
				sum_nu_kappa += b[nu]*b[kappa]*moment6(times);
			}
			sum_kappa += b[kappa]*(moment4(j + i - kappa) + moment4(j - kappa));
		}
		Mi(j, 0) = (n - q - i) * (sum_nu_kappa + a_zero*sum_kappa + a_zero*a_zero*coeff[EG2c]) / (n - q); // half of 3.14
	}

	// k > 0, j = 0
	for (int k = 1; k <= q; k++) {
		scalar sum_nu_kappa = 0;
		scalar sum_kappa = 0;
		for (int kappa = 0; kappa <= q; kappa++) {
			for (int nu = 0; nu <= q; nu++) {
				int time[3];
				time[0] = i -kappa;
				time[1] = -nu;
				time[2] = i - k;
				sum_nu_kappa += b[nu]*b[kappa]*moment6(time);
			}
			sum_kappa += b[kappa]*(moment4(kappa - k) + moment4(kappa + i - k));
		}
		Mi(0, k) = (n - q - i) * (sum_nu_kappa + a_zero*sum_kappa + a_zero*a_zero*coeff[EG2c]) / (n - q); // half of 3.14
	}

	// j > 0, k > 0
	for (int j = 1; j <= q; j++) {
		for (int k = 1; k <= q; k++) {
			scalar sum_nu_kappa = 0;
			scalar sum_nu = 0;
			for (int nu = 0; nu <= q; nu++) {
				for (int kappa = 0; kappa <= q; kappa++) {
					int time[4];

					// in sorensen r occurs in every time but since we only need the differences we can remove it
					time[0] = - nu;
					time[1] = i - kappa;
					time[2] = - j;
					time[3] = i - k;

					sum_nu_kappa += b[nu]*b[kappa]*moment8(time);
				} // kappa

				int timeA[3];
				// we removed r from the times
				timeA[0] = i - nu;
				timeA[1] = - j;
				timeA[2] = i - k;
				int timeB[3];
				timeB[0] = -nu;
				timeB[1] = -j;
				timeB[2] = i - k;
                                
				sum_nu += b[nu]*(moment6(timeA) + moment6(timeB));
			} // nu

			Mi(j, k) = (n - q - i) * (sum_nu_kappa + a_zero*sum_nu + a_zero*a_zero*moment4(k - j - i)) / (n - q); // half of 3.14
		} // k
	} // j
}

void computeMi_R(
		const int *iParam,
		const int *qParam,            // q
		const int *nParam,            // n
		const double *hParam,         // r
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		const double b[],              // 
		const double *a_zeroParam,    //
		double Mi_result[]) {

	scalar *psiexp = computePsiExp(nParam, psiParam);
	//scalar *coeff = computeCoefficients(hParam, thetaParam, psiParam);
	#include "coefficientsAsArray.part.c"

	const int q = *qParam;
	scalar *Mi = malloc(sizeof(scalar)*(q + 1)*(q + 1));
	computeMi(*iParam, *qParam, *nParam, b, *a_zeroParam, psiexp, coeff, Mi);

	for (int jk = 0; jk < (q + 1)*(q + 1); jk++) {
		Mi_result[jk] = (double)Mi[jk];
	}

	free(Mi);
	free(psiexp);
	free(coeff);
}

void computeMloop(
		const int *qParam,            // q
		const int *nParam,            // n
		const double *hParam,         // t
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		const double b[],              // 
		const double *a_zeroParam,    //
		double M_result[]              // the resulting (q + 1) square matrix
		) {
	// assign parameters to variables
	const int q = *qParam;
	const int n = *nParam;
	if (n < q + 1) {
		printf("computeM: parameter n = %d is smaller than (q + 1 = %d)\n", n, q + 1);
		return;
	}

	// precompute exponentials and coefficients
	scalar *psiexp = computePsiExp(nParam, psiParam);
	#include "coefficientsAsArray.part.c"

	scalar *Macc = malloc(sizeof(scalar)*(q + 1)*(q + 1));
	scalar *Mi = malloc(sizeof(scalar)*(q + 1)*(q + 1));

	printf("Computing M[0]... ");
	fflush(stdout);
	computeMi(0, *qParam, *nParam, b, *a_zeroParam, psiexp, coeff, Macc);
	printf("done.\n");
	for (int i = 1; i <= n - q - 1; i++) {
		printf("Computing M[%d]... ", i);
		fflush(stdout);
		computeMi(i, *qParam, *nParam, b, *a_zeroParam, psiexp, coeff, Mi);
		printf("done.\n");
		for (int j = 0; j <= q; j++) {
			for (int k = 0; k <= q; k++) {
				Macc[j*(q + 1) + k] += Mi[j*(q + 1) + k] + Mi[k*(q + 1) + j];
			}
		}
	}

	for (int jk = 0; jk < (q + 1)*(q + 1); jk++) {
		M_result[jk] = (double)Macc[jk];
	}

	free(Macc);
	free(Mi);
	free(psiexp);
	free(coeff);
}

void computeMklu(
		const int *dParam,            // d
		const int *nParam,            // n
		const double *hParam,         // t
		const double thetaParam[4],   // the parameters, in order: beta, eta, tau, phi
		const double psiParam[4],     // psi1, psi2, psi3, psi4
		double M_result[]              // the resulting (q + 1) square matrix
		) {
	// assign parameters to variables
	const int d = *dParam;
	const int n = *nParam;
	const int nd = 2*n;

	// precompute exponentials and coefficients
	scalar *psiexp = computePsiExp(&nd, psiParam);
	#include "coefficientsAsArray.part.c"
	
	{
		        scalar sum_i = 0;
                for (int i = 1; i <= n; i++) {
                    sum_i += moment4(i);
				}

				M_result[0] = (double)(moment4(0) + 2*sum_i - (2*n + 1)*coeff[EG2c]*coeff[EG2c]);
            
        }
	
	for (int k = 0; k <= d; k++) {
		        scalar sum_i = 0;
                for (int i = 1; i <= n; i++) {
                    int times[3];
                    times[0] = 1;
                    times[1] = 1 + i;
                    times[2] = 1 + k;
                    sum_i += moment6(times);
					times[0] = 1;
                    times[1] = 1 + i;
                    times[2] = 1 + i + k;
                    sum_i += moment6(times);
                }
				
                int times[3];
                times[0] = 1;
                times[1] = 1 + k;
                times[2] = 1;

				M_result[k + 1] = M_result[(k + 1)*(d + 2)] = (double)(moment6(times) + sum_i - (2*n + 1)*moment4(k)*coeff[EG2c]);
            
        }

	
        for (int k = 0; k <= d; k++) {
            for (int l = 0; l <= d; l++) {
                scalar sum_i = 0;
                for (int i = 1; i <= n; i++) {
                    int times[4];
                    times[0] = 1;
                    times[1] = 1 + i;
                    times[2] = 1 + k;
                    times[3] = 1 + i + l;
                    sum_i += moment8(times);
					times[0] = 1;
                    times[1] = 1 + i;
                    times[2] = 1 + l;
                    times[3] = 1 + i + k;
                    sum_i += moment8(times);
                }
				
                int times[4];
                times[0] = 1;
                times[1] = 1;
                times[2] = 1 + k;
                times[3] = 1 + l;

                M_result[(k + 1)*(d + 2) + l + 1] = (double)(moment8(times) + sum_i - (2*n + 1)*moment4(k)*moment4(l));
            }
        }

	free(psiexp);
	free(coeff);
}


