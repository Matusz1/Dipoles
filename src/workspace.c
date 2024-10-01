#include <string.h>

#include <lapacke.h>

#include "workspace.h"
#include "fock.h"
#include "state.h"
#include "util.h"

Workspace workspace_create_from_basis(FockBasis* basis) {
        Workspace w = {
                .basis = *basis,
                .npot = 0,
                .V = { (void*) 0 },
                .pdat = { (void*) 0 },
        };
        basis->size = 0;
        basis->cap = 0;
        basis->states = (void*) 0;

        const uint s = w.basis.size;
        w.H0 = util_malloc(sizeof(*w.H0)*s);
        for (uint i = 0; i < s; ++i)
                w.H0[i] = fock_energy(&w.basis.states[i]);

        return w;
}

void workspace_delete(Workspace *w) {
        fock_basis_free(&w->basis);
        util_free(w->H0);
        for (uint i = 0; i < w->npot; ++i) {
                state_potential_data_free(w->pdat[i]);
                util_free(w->V[i]);
        }
}

void workspace_add_potential(Workspace* w, const char* fname) {
        util_error(w->npot == 10, "Cannot add more potentials.\n");
        w->pdat[w->npot] = state_potential_alloc_from_file(fname);

        const uint s = w->basis.size;
        w->V[w->npot] = util_malloc(sizeof(*w->V[w->npot])*s*s);
        // TODO: Computig the potental should be done once and for all,
        // not one by one after every call to workspace_add_potential(...)
        for (uint i = 0; i < s; ++i) {
                for (uint j = i; j < s; ++j) {
                        double v = fock_potential(&w->basis.states[i], &w->basis.states[j], w->pdat[w->npot]);
                        w->V[w->npot][i*s+j] = v;
                        w->V[w->npot][j*s+i] = v;
                }
        }
        
        ++w->npot;
}

double* __workspace_alloc_hamiltonian(Workspace* w, double* g) {
        util_error(w->basis.size == 0, "Size of basis is equal to zero.");

        const uint s = w->basis.size;
	double* H = util_malloc(sizeof(*H)*s*s);
        memset(H, 0, sizeof(*H)*s*s);

	for (uint i = 0; i < s; ++i) {
		double t = w->H0[i];
                H[i*s + i] = t;
		for (uint j = i; j < s; ++j) {
                        double v = 0.0;
                        for (uint k = 0; k < w->npot; ++k)
                                v += g[k] * w->V[k][i*s + j];
			H[i*s + j] += v;
		}
	}
        return H;
}

DEigenSol workspace_compute(Workspace* w, double* g) {
        /*
        uint max = 0;
        for (uint i = 0; i < s; ++i)
                for (uint j = 0; j < basis.states[i].size; ++j)
                        max = MAX(basis.states[i].states[j], max);
        */
        util_error(w->basis.size == 0, "Size of basis is equal to zero.");

        // TODO: Add fucntion workspace_alloc_hamiltonian
        const uint s = w->basis.size;
	DEigenSol ret = {
		.size = s,
		.U = __workspace_alloc_hamiltonian(w, g),
		.val = util_malloc(sizeof(*ret.val)*s)
	};

	LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', s, ret.U, s, ret.val);
	return ret;
}

DEigenSol workspace_compute_Hcm(Workspace* w, double* g, double* Hcm) {
        util_error(w->basis.size == 0, "Size of basis is equal to zero.");

        // TODO: Add fucntion workspace_alloc_hamiltonian
        const uint s = w->basis.size;
	DEigenSol ret = {
		.size = s,
		.U = __workspace_alloc_hamiltonian(w, g),
		.val = util_malloc(sizeof(*ret.val)*s)
	};

	for (uint i = 0; i < s; ++i)
		for (uint j = i; j < s; ++j)
                        ret.U[i*s + j] += Hcm[i*s + j];
	LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', s, ret.U, s, ret.val);
	return ret;
}

void deigen_sol_free(DEigenSol *es) {
        free(es->U);
        free(es->val);
        es->size = 0;
}
