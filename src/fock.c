#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <lapacke.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_heapsort.h>

#include "fock.h"
#include "state.h"
#include "util.h"

static inline bool is_even(int n) {
        return (n % 2) == 0;
}

static inline bool is_odd(int n) {
        return !is_even(n);
}

Fock fock_init(uint size) {
	util_error(size > CONFIG_FOCK_CAPACITY, "fock_init, size too big\n");
        Fock fock = {
                .size = size,
                .states = { 0 },
        };
	return fock;
}

bool fock_equal(const Fock* s1, const Fock* s2) {
        return (s1->size == s2->size)
                && !memcmp(s1->states, s2->states, s1->size*sizeof(*s1->states));
}



/* Functions for dealing with fock basis */
/* ===================================== */

void fock_basis_add(FockBasis* set, const Fock* fock) {
	++set->size;
        if (set->cap < set->size) {
                set->cap = MAX(set->size, 2*set->cap);
	        set->states = util_realloc(set->states, sizeof(*set->states) * set->cap);
        }
	memcpy(set->states + (set->size - 1), fock, sizeof(*fock));
}

void __fock_basis_alloc_recursive(FockBasis* set, Fock* fock, uint npart, double e_cutoff, bool accept_fun(Fock*, void*), void* args) {
        if (npart == 0) {
                if (accept_fun(fock, args))
                        fock_basis_add(set, fock);
                return;
        }

        for (uint i = 0; i <= fock->states[npart]; ++i) {
                fock->states[npart-1] = i;
                double next_cutoff = e_cutoff - state_energy(i);
                if (next_cutoff >= 0)
                        __fock_basis_alloc_recursive(set, fock, npart-1, e_cutoff, accept_fun, args);
        }
        // This is probably not needed, it used to be required
        //fock->states[npart-1] = 0; /* Reset to the least energy state */
}


FockBasis fock_basis_alloc_special(uint npart, double e_cutoff, bool accept_fun(Fock*, void*), void* args) {
        FockBasis basis = { 0 };
        Fock fock = {
                .size = npart,
                .states = { 0 },
        };

        uint max = 0;
        while (state_energy(max) < e_cutoff)
                ++max;

        // TODO: There might be some additional
        // cutoff on the states, something like this should be added.
        //max = MIN(max, HERMITE_MAX_ORDER+1);

        for (uint i = 0; i < max; ++i) {
                fock.states[npart-1] = i;
                double next_cutoff = e_cutoff - state_energy(i);
                if (next_cutoff >= 0)
                        __fock_basis_alloc_recursive(&basis, &fock, npart-1, next_cutoff, accept_fun, args);
        }
        return basis;
}

static bool __fock_true(Fock* fock, void* args) {
        return true;
}

FockBasis fock_basis_alloc(uint npart, double e_cutoff) {
        return fock_basis_alloc_special(npart, e_cutoff, __fock_true, (void*) 0);
}

double fock_energy(const Fock* fock) {
	double E = 0.0;
	for (uint i = 0; i < fock->size; ++i)
		E += state_energy(fock->states[i]);
	return E;
}

void fock_basis_free(FockBasis* set) {
	free(set->states);
	set->size = 0;
        set->cap = 0;
}

uint fock_operator_count(const Fock* fock, uint nr) {
	uint count = 0;
	for (uint i = 0; i < fock->size; ++i)
		if (fock->states[i] == nr)
			++count;
	return count;
}

double fock_operator_create(Fock* fock, uint nr) {
	util_error(fock->size == CONFIG_FOCK_CAPACITY, "Cannot add more particles to fock state.\n");
	uint idx = fock->size;
        for (; idx > 0; --idx)
                if (nr < fock->states[idx-1])
                        fock->states[idx] = fock->states[idx-1];
                else
                        break;
	fock->states[idx] = nr;
	++fock->size;
	return sqrt(fock_operator_count(fock, nr));
}

double fock_operator_annihilate(Fock* fock, uint nr) {
	uint idx = 0;
	uint count = 0;
	for (uint i = 0; i < fock->size; ++i) {
		if (nr == fock->states[i]) {
			++count;
			idx = i;
		}
	}
	/* Special case when we cannot annichilate any particles,
	 * state just becomes 0 */
	if (count == 0) {
		fock->size = 0;
		return 0.0;
	}
	memmove(fock->states+idx, fock->states+idx+1, sizeof(*fock->states)*(fock->size-idx-1));
	--fock->size;
	return sqrt((double)count);
}

double fock_potential(const Fock* left, const Fock* right, PotentialData* dat) {
        // TODO: Before we even start computing the matrix elements,
        // we can check if the states can possibly be nonzero,
        // that is if they differ no more than by two single part. states.
        
	Fock a1;
	Fock a2;
	Fock a3;
	Fock a4;

	double ret_v = 0.0;

	/* This is a bit tricky because we have to iterate over every UNIQUE state */
	uint p1 = right->states[0];
	uint p2 = right->states[0];
	uint p3 = left->states[0];
	uint p4 = left->states[0];

	for (uint i = 0; i < right->size; ++i) {
		if (i != 0 && p1 == right->states[i])
			continue;
                a1 = *right;
		double v1 = fock_operator_annihilate(&a1, right->states[i]);
		p1 = right->states[i];
		for (uint j = 0; j < a1.size; ++j) {
			if (j != 0 && p2 == a1.states[j])
				continue;
                        a2 = a1;
			double v2 = fock_operator_annihilate(&a2, a1.states[j]);
			p2 = a1.states[j];
			for (uint q = 0; q < left->size; ++q) {
				if (q != 0 && p3 == left->states[q])
					continue;
                                a3 = a2;
				double v3 = fock_operator_create(&a3, left->states[q]);
				p3 = left->states[q];
				for (uint p = 0; p < left->size; ++p) {
					if (p != 0 && p4 == left->states[p])
						continue;
                                        a4 = a3;
					double v4 = fock_operator_create(&a4, left->states[p]);
					p4 = left->states[p];
					if (fock_equal(&a4, left)) {
                                                double v = state_potential(p3, p4, p1, p2, dat);
						ret_v += v1*v2*v3*v4*v;
                                        }

				}
			}
		}
	}

        /* Coefficient comes from 1/2 before sum in interaction term */
        return 0.5 * ret_v;
}

static void __dec_to_binarr(int* res, uint n, int dim) {
	int pos = dim - 1;
	for (uint i = 0; i < dim+1; ++i)
		res[i] = 0;

	while (n) {
		res[pos] = n & 1;
		res[dim] += res[pos];
                n = n >> 1;
		//n = n / 2;
		--pos;
	}
}

#define MINUS_ONE_POW(n) (((n)%2)==0 ? 1 : -1)

double __permanent(double* A, int n) {
	static int chi[CONFIG_FOCK_CAPACITY+1];
        uint C = (1 << n);
	//const double C = pow(2, n); 
	double sum = 0;
	double rowsumprod, rowsum;

	for (uint k = 1; k < C; ++k) {
		rowsumprod = 1;
		__dec_to_binarr(chi, k, n);

		for (int m = 0; m < n; ++m) {
			rowsum = 0;
			for (int p = 0; p < n; ++p)
				rowsum += chi[p] * A[m * n + p];
			rowsumprod *= rowsum;    
		}        

		sum += MINUS_ONE_POW(n-chi[n]) * rowsumprod;
		//sum += pow(-1.0, n - chi[n]) * rowsumprod;
	}    

	return sum;
}

// TODO: Make some CONFIG constant that will allow
// making taking more dimensions

double fock_compute(Fock* fock, double* x) {
        const uint s = fock->size;
	double A[s*s];
	uint ij = 0;

	for (uint i = 0; i < s; ++i)
		for (uint j = 0; j < s; ++j, ++ij)
			A[ij] = state_val(fock->states[j], x[i]);


        double norm_fact = gsl_sf_fact(s);
        /* Count how many times the state repeats */
        for (uint i = 0; i < s;) {
                const uint last = i;
                while (i < s && fock->states[last] == fock->states[i])
                        ++i;
                norm_fact *= gsl_sf_fact(i - last);
        }

	return __permanent(A, s) / sqrt(norm_fact);
}

double fock_basis_compute(FockBasis* basis, double* a, double* x) {
        double sum = 0.0;
        for (uint i = 0; i < basis->size; ++i)
                sum += a[i] * fock_compute(&basis->states[i], x);
        return sum;
}

/* === Printing === */
/* ================ */

void fock_fprint(FILE* f, const Fock* fock) {
	for (uint i = 0; i < fock->size;) {
		fprintf(f, "%u", fock->states[i]);
		if (++i != fock->size)
			putc(' ', f);
	}
}

void fock_basis_fprint(FILE* f, const FockBasis set) {
	for (uint i = 0; i < set.size; ++i) {
		fprintf(f, "%6u: ", i);
		fock_fprint(f, &set.states[i]);
		putc('\n', f);
	}
}

void fock_basis_fprint_weighted(FILE* f, const FockBasis set, const double* coeff, uint n) {
	size_t idxs[set.size];
	//gsl_heapsort_index(idxs, coeff, set.size, sizeof(*coeff), __fock_square_compare);
	fprintf(f, " |ci|^2         state\n");
	for (uint i = 0; i < MIN(set.size, n); ++i) {
		fprintf(f, "%lf : ", coeff[idxs[i]]*coeff[idxs[i]]);
		fock_fprint(f, &set.states[idxs[i]]);
		putc('\n', f);
	}
}
