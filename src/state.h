#ifndef _STATE_H_
#define _STATE_H_

#include "core.h"

double state_val(uint n, double x);

static inline double state_energy(uint n) {
        return 0.5 + n;
}

typedef struct PotentialData PotentialData;

PotentialData* state_potential_alloc_from_file(const char* fname);
void state_potential_data_free(PotentialData* pdat);

double state_potential(int m1, int m2, int m3, int m4, PotentialData* dat);

typedef struct Fock Fock;

/* Function used as argument to 'fock_basis_alloc_special',
 * it decides whether the state has given parity.
 * 'even' point to boolean with true if we want to check if 
 * the fock state is even. */
bool fock_hermite_parity(Fock* fock, void* even);


/* Two particles computations */
void example_center_of_mass(double gdd);
/* This uses vdW units for Dy */
void example_center_of_mass_over_a(double gdd);
void example_center_of_mass_print_wf(double g, double gdd);
void example_convergence_over_a(double oa);


#endif // _STATE_H_
