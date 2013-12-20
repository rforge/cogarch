#include <math.h>

#include "common.h"

static inline scalar Ceiling(const scalar x) {
#ifdef LD
	return ceill(x);
#else
	return ceil(x);
#endif
}

static inline scalar Floor(const scalar x) {
#ifdef LD
	return floorl(x);
#else
	return floor(x);
#endif
}
	

static inline int Max(const int i, const int j) {
	if (i >= j)
		return i;
	else
		return j;
}

static inline int Max3(const int i, const int j, const int k) {
    return Max(i, Max(j, k));
}

static inline int Min(const int i, const int j) {
	if (i <= j)
		return i;
	else
		return j;
}

static inline int Min3(const int i, const int j, const int k) {
    return Min(i, Min(j, k));
}

static inline scalar Exponential(const scalar x) {
#ifdef LD
	return expl(x);
#else
	return exp(x);
#endif
}

static inline scalar Power(const scalar b, const scalar e) {
#ifdef LD
	return powl(b, e);
#else
	return pow(b, e);
#endif
}

static inline scalar Cosh(const scalar x) {
#ifdef LD
	return coshl(x);
#else
	return cosh(x);
#endif
}

static inline scalar Sinh(const scalar x) {
#ifdef LD
	return sinhl(x);
#else
	return sinh(x);
#endif
}

// Compute powers with fixed integer exponents as multiplications. 
// The number of multiplications should be minimal.
static inline scalar Power2(scalar x) { return x*x; }			// 1 mult
static inline scalar Power3(scalar x) { return Power2(x)*x; }		// 2 mult
static inline scalar Power4(scalar x) { return Power2(Power2(x)); }	// 2 mult
static inline scalar Power5(scalar x) { return Power4(x)*x; }		// 3 mult
static inline scalar Power6(scalar x) { return Power3(Power2(x)); }	// 3 mult
static inline scalar Power7(scalar x) { return Power6(x)*x; }		// 4 mult
static inline scalar Power8(scalar x) { return Power4(Power2(x)); }	// 3 mult
static inline scalar Power9(scalar x) { return Power8(x)*x; }		// 4 mult
static inline scalar Power10(scalar x) { return Power5(Power2(x)); }	// 4 mult
static inline scalar Power11(scalar x) { return Power10(x)*x; }		// 5 mult
static inline scalar Power12(scalar x) { return Power6(Power2(x)); }	// 4 mult
static inline scalar Power13(scalar x) { return Power12(x)*x; }		// 5 mult
static inline scalar Power14(scalar x) { return Power7(Power2(x)); }	// 5 mult
static inline scalar Power15(scalar x) { return Power14(x)*x; }		// 6 mult
static inline scalar Power16(scalar x) { return Power8(Power2(x)); }	// 4 mult


