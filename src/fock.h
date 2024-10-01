#ifndef _FOCK_H_
#define _FOCK_H_

#include <stdio.h>

#include "core.h"

/* Max number of particles in fock state */
#define CONFIG_FOCK_CAPACITY 6

typedef struct Fock {
	uint states[CONFIG_FOCK_CAPACITY]; /* Single particles states in sorted order 0, ..., n */
	uint size; /* Number of particles */
} Fock;

Fock fock_init(uint size);

bool fock_equal(const Fock* s1, const Fock* s2);

/* === Way of handling sets of fock states === */
/* =========================================== */

/* This structers allows indexing Fock states,
 * f.e. constructing linear combinations for same basis */

typedef struct {
	Fock* states;
	uint size;
        uint cap;
} FockBasis;

/* Generate basis according to energy cutoff and function
 * that can either accept or reject given state. 
 * For example if you want states of say given parity
 * 'accept_fun' decides whether state has given parity, 
 * using 'args' you can pass arguments to 'accept_fun' */
FockBasis fock_basis_alloc_special(uint npart, double e_cutoff, bool accept_fun(Fock*, void*), void* args);

/* Construct basis with just the energy cutoff requirement */
FockBasis fock_basis_alloc(uint npart, double e_cutoff);

double fock_energy(const Fock* fock);

void fock_basis_add(FockBasis* basis, const Fock* fock);
void fock_basis_free(FockBasis* basis);

uint fock_operator_count(const Fock* fock, uint nr);
double fock_operator_create(Fock* fock, uint nr);
double fock_operator_annihilate(Fock* fock, uint nr);

/* Way of calculating interaction potential between two states */

typedef struct PotentialData PotentialData;
double fock_potential(const Fock* left, const Fock* right, PotentialData* pot);

double fock_compute(Fock* fock, double* x);
double fock_basis_compute(FockBasis* basis, double* a, double* x);


/* === Printing === */
/* ================ */

void fock_fprint(FILE* f, const Fock* fock);
static inline void fock_print(const Fock* fock) {
	fock_fprint(stdout, fock);
}

void fock_basis_fprint(FILE* f, const FockBasis set);
static inline void fock_basis_print(const FockBasis set) {
	fock_basis_fprint(stdout, set);
}

void fock_basis_fprint_weighted(FILE* f, const FockBasis set, const double* coeff, uint n);
static inline void fock_basis_print_weighted(const FockBasis set, const double* coeff, uint n) {
	fock_basis_fprint_weighted(stdout, set, coeff, n);
}

/*
void fock_coeffstateset_fprint(FILE* f, const FockState set, uint n);
static inline void fock_coeffstateset_print(const FockState set, uint n) {
	fock_coeffstateset_fprint(stdout, set, n);
}
*/

#endif // _FOCK_H_
