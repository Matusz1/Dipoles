#include "hermite.h"
#include <math.h>
#include <string.h>

#include "util.h"
#include "polynomial.h"

/* We are precalculate the hermite coefficients using
 * modified recurence formula (modified, meaning we keep the
 * polynomials normalized) */

static uint _hermite_max_deg = 0;
static double* _hermite_coeff = (void*) 0;

void _hermite_precalc(uint max_deg) {
        util_error(max_deg > HERMITE_MAX_ORDER, "Max degree of Hermite polynomial is %u.\n", HERMITE_MAX_ORDER);
        max_deg = MAX(1, max_deg);
        if (_hermite_max_deg >= max_deg)
                return;

        // If we want hermite polynomial of degree n available,
        // we need n+1 polynomials (0,1,2,...,n)
        uint s = max_deg + 1;
        _hermite_coeff = util_realloc(_hermite_coeff, sizeof(*_hermite_coeff)*s*s);
        memset(_hermite_coeff, 0, sizeof(*_hermite_coeff)*s*s);

        // Coeffients for 0th and 1st order Hermite
        _hermite_coeff[0]   = 1.0 / sqrt(sqrt(M_PI));
        _hermite_coeff[s]   = 0.0;
        _hermite_coeff[s+1] = sqrt(2.0/sqrt(M_PI));

        // Determine the rest using recurence
        for (int n = 1; n < s-1; ++n) {
                double fact = sqrt(2.0/(n+1));
                for (uint j = 0; j <= n; ++j)
                        _hermite_coeff[s*(n+1)+j+1] = fact*_hermite_coeff[s*n+j];
                fact = -sqrt((double) n/(n+1));
                for (uint j = 0; j < n; ++j)
                        _hermite_coeff[s*(n+1)+j] += fact*_hermite_coeff[s*(n-1)+j];
        }

        // Alternative way using explicitly the coefficents
        /*
        for (int n = 0; n <= max_deg; ++n) {
                double fact = 1.0 / sqrt(sqrt(M_PI));
                for (int i = 1; i <= n; ++i)
                        fact *= sqrt(i/2.0);
                for (int i = 0; i <= n/2; ++i) {
                        double a = (i%2 == 0) ? 1 : -1;
                        for (int j = 1; j <= n-2*i; ++j)
                                a *= 2.0 / j;
                        for (int j = 1; j <= i; ++j)
                                a /= j;
                        _hermite_coeff[n*s+(n-2*i)] = a*fact;
                }
        }
        */
        _hermite_max_deg = max_deg;
}

Polynomial hermite_pol(uint n) {
        _hermite_precalc(n);
        Polynomial p = {
                .cap = -1,
                .deg = n,
                .a = &_hermite_coeff[n*(_hermite_max_deg+1)],
        };
        return p;
}

double hermite_val(uint n, double x) {
        Polynomial p = hermite_pol(n);
        return polynomial_val(&p, x);
}
