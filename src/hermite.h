#ifndef _HERMITE_H_
#define _HERMITE_H_

#include "core.h"
#include "polynomial.h"

#define HERMITE_MAX_ORDER 50

/* IMPORTANT: In this code, hermite polynomials are normalized
 * with respect to usual integral with weight exp(-x*x). */

Polynomial hermite_pol(uint n);
double hermite_val(uint n, double x);

#endif // !_HERMITE_H_
