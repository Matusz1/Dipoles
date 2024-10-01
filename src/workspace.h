#ifndef _WORKSPACE_H_
#define _WORKSPACE_H_

#include "core.h"
#include "fock.h"

/* The computation workspace is used to hold matrices
 * of noninteracting Hamiltonian, different interaction
 * potentials and data that might have been precomputed
 * and lying in outer file, that helps determine those
 * matrices (PotentialData). By having this structure we can 
 * diagonalize Hamiltonian without having to recalculate 
 * the whole matrix again, just use different coupling constant. */

typedef struct Workspace {
        FockBasis basis;

        /* Diagonal of the noninteracting Hamiltonian */
        double* H0;

        uint npot;
        double* V[10];
        PotentialData* pdat[10];
} Workspace;

/* This function takes the ownership of a basis, resets the pointer */
Workspace workspace_create_from_basis(FockBasis* basis);
void workspace_delete(Workspace *w);
// TODO: Function return Hamiltonian

void workspace_add_potential(Workspace* w, const char* fname);

typedef struct DEigenSol {
        double* U;
        double* val;
        uint size;
} DEigenSol;

void deigen_sol_free(DEigenSol *es);

DEigenSol workspace_compute(Workspace* w, double* g);

// TODO: Better way to use a workspace, we need workspace_alloc_hamiltonian()
DEigenSol workspace_compute_Hcm(Workspace* w, double* g, double* Hcm);


#endif // !_WORKSPACE_H_
