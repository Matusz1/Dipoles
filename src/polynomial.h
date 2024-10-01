#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "core.h"

typedef struct Polynomial {
        uint deg;  // Degree of a polynomial
        int cap;   // Capacity of memory for a, for maximum degree,
                   // if negative a is unchangable
        double *a;
} Polynomial;

double polynomial_val(Polynomial *p, double x);

typedef struct Pol2 {
        uint deg;
        int cap;
        double *a;
} Pol2;

void pol2_free(Pol2* p);
Pol2 pol2_from_p(Polynomial* p);

#endif // _POLYNOMIAL_H_
