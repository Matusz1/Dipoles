#include "polynomial.h"
#include "util.h"

double binomial(uint n, uint k) {
        if (n == k || k == 0)
                return 1.0;
        double v = n;
        for (uint i = 1; i < k; ++i)
                v *= (double)(n-i)/i;
        return v;
}

double polynomial_val(Polynomial *p, double x) {
        double v = p->a[p->deg];
        for (uint i = p->deg; i > 0; --i)
                v = p->a[i-1] + x*v;

        return v;
}

void pol2_free(Pol2* p) {
        p->deg = 0;
        if (p->cap >= 0)
                util_free(p->a);
}

Pol2 pol2_from_p(Polynomial* p) {
        Pol2 ret = {
                .deg = p->deg,
                .cap = p->deg,
        };
        uint s = p->deg+1;
        ret.a = util_malloc(sizeof(*ret.a)*s*s);

        for (uint i = 0; i <= p->deg; ++i)
                for (uint j = 0; j <= p->deg; ++j)
                        ret.a[i*s+j] = p->a[i] * binomial(i, j);
        return ret;
}
