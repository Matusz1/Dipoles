#include <math.h>
#include <stdio.h>

#include "hermite.h"
#include "polynomial.h"

void test_hermite() {
        const uint npts = 1000;
        for (uint n = 0; n <= HERMITE_MAX_ORDER; ++n) {
                const double lim = 2.0*sqrt(n+10);
                const double h = 2*lim / npts;

                double sum = 0.0;
                for (uint i = 1; i < npts-1; ++i) {
                        const double x = -lim + i*h;
                        const double e = exp(-x*x);
                        const double v = hermite_val(n, x);
                        sum += v*v*e;
                }
                double v = hermite_val(n, -lim);
                sum += 0.5*v*v*exp(-lim*lim);
                v = hermite_val(n, lim);
                sum += 0.5*v*v*exp(-lim*lim);
                sum *= h;
                printf("%2u   %.16lf\n", n, sum);
        }
}

void test_pol2() {
        // 1 + 2*x + 3*x*x
        double a[3] = {1, 2, 3};
        Polynomial p = {
                .deg = 2,
                .cap = -1,
                .a = a,
        };
        Pol2 p2 = pol2_from_p(&p);
        for (uint i = 0; i < 3; ++i) {
                for (uint j = 0; j < 3; ++j) {
                        printf("%6.2lf ", p2.a[j*3 + i]);
                }
                putchar('\n');
        }
        pol2_free(&p2);
}
