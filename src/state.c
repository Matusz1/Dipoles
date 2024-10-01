#include <tgmath.h>
#include <string.h>

#include "fock.h"
#include "util.h"
#include "hermite.h"
#include "workspace.h"

struct PotElem {
        u8 idx[4];
        double v;
};

/* 'Outer world' only knows that such structure exists,
 * precise implementation is dependent on the (here) oscillator states */
typedef struct PotentialData {
        uint size;
        struct PotElem* data;
} PotentialData;

double state_val(uint n, double x) {
        return  exp(-0.5*x*x)*hermite_val(n, x);
}

PotentialData* state_potential_alloc_from_file(const char* fname) {
        FILE* f = util_fopen(fname, "rb");

        u32 s;
        util_fread(&s, sizeof(s), 1, f);

        PotentialData* ret = util_malloc(sizeof(*ret));
        ret->data = util_malloc(sizeof(*ret->data)*s);

        for(uint i = 0; i < s; ++i) {
                util_fread(ret->data[i].idx, sizeof(ret->data->idx), 1, f);
                util_fread(&ret->data[i].v, sizeof(ret->data->v), 1, f);
        }

        util_fclose(f);
        ret->size = s;

        return ret;
}

void state_potential_data_free(PotentialData* pdat) {
        util_free(pdat->data);
        util_free(pdat);
}

static int __pdata_comp(const void* lhs, const void* rhs) {
        const u8* p1 = ((const struct PotElem*) lhs)->idx;
        const u8* p2 = ((const struct PotElem*) rhs)->idx;

        // FIXME: Maybe there is a better solution?
        if (p1[3] < p2[3])
                return -1;
        else if (p1[3] > p2[3])
                return  1;

        if (p1[2] < p2[2])
                return -1;
        else if (p1[2] > p2[2])
                return  1;

        if (p1[1] < p2[1])
                return -1;
        else if (p1[1] > p2[1])
                return  1;

        if (p1[0] < p2[0])
                return -1;
        else if (p1[0] > p2[0])
                return  1;

        return 0;
}

static inline void __swap_u8(u8* l, u8* r) {
        u8 t = *r;
        *r = *l;
        *l = t;
}

double state_potential(int m1, int m2, int m3, int m4, PotentialData* dat) {
        if ((m1+m2+m3+m4) % 2)
                return 0.0;

        /* First make sure the second number in each pair (1,2), (3,4) is bigger */
        u8 key[4] = { MIN(m1,m2), MAX(m1,m2), MIN(m3,m4), MAX(m3,m4) };
        /* Make sure the biggest number is at the last position */
        if (key[3] < key[1]) {
                __swap_u8(key, key+2);
                __swap_u8(key+1, key+3);
        }
        /* Special case when (2) == (4), then it must be (1) <= (3) */
        if (key[1] == key[3] && key[0] > key[2])
                __swap_u8(key, key+2);

        struct PotElem* e = bsearch(key, dat->data, dat->size, sizeof(*dat->data), __pdata_comp);
        if (e)
                return e->v;

        util_error(true, "Failed to fill potential!\n"
                   "Serching failed for: (%2d|%2d|%2d|%2d)\n",
                   key[0], key[1], key[2], key[3]);
}

bool fock_hermite_parity(Fock* fock, void* even) {
        uint tot = 0;
        for (uint i = 0; i < fock->size; ++i)
                tot += fock->states[i];
        return (*((bool*)even) == !(tot % 2));
}

/* Lots of two particles computations with some weird functions */

/* Here we are computing center of mass hamiltonian for two particles,
 * probably not the most elegant solution, but it works! */

static double __Hcm_int_dx(uint m, uint n) {
        if (m == n+1)
                return -sqrt(n+1);
        if (m+1 == n)
                return sqrt(n);
        return 0.0;
}

static double __Hcm_int_x(uint m, uint n) {
        if (m == n+1)
                return sqrt(n+1);
        if (m+1 == n)
                return sqrt(n);
        return 0.0;
}

static double __Hcm_int(uint mx1, uint my1, uint my2, uint mx2) {
        double sum = 0.0;
        sum -= __Hcm_int_dx(mx1, mx2) * __Hcm_int_dx(my1, my2);
        sum += __Hcm_int_x(mx1, mx2) * __Hcm_int_x(my1, my2);
        // TODO: Check if this really should be 8, definitely Ecm = 0.5
        return sum / 8;
}

static double __Hcm_mval(Fock* fl, Fock* fr) {
        const uint* ul = fl->states;
        const uint* ur = fr->states;
        double sum = 0.0;
        if (ur[0] == ur[1]) {
                if (ul[0] == ul[1]) {
                        return 2*__Hcm_int(ul[0], ul[1], ur[0], ur[1]);
                }
                sum += __Hcm_int(ul[0], ul[1], ur[0], ur[1]);
                sum += __Hcm_int(ul[1], ul[0], ur[0], ur[1]);
                return sqrt(2)*sum;
        }
        if (ul[0] == ul[1]) {
                sum += __Hcm_int(ul[0], ul[1], ur[0], ur[1]);
                sum += __Hcm_int(ul[0], ul[1], ur[1], ur[0]);
                return sqrt(2)*sum;
        }
        sum += __Hcm_int(ul[0], ul[1], ur[0], ur[1]);
        sum += __Hcm_int(ul[0], ul[1], ur[1], ur[0]);
        sum += __Hcm_int(ul[1], ul[0], ur[0], ur[1]);
        sum += __Hcm_int(ul[1], ul[0], ur[1], ur[0]);

        return sum;
}

static void __Hcm_square(double* Hcm, uint s) {
        for (uint i = 0; i < s; ++i)
                Hcm[i*s + i] -= 0.5;

        double* Hnew = util_malloc(sizeof(*Hnew)*s*s);
        memset(Hnew, 0, sizeof(*Hnew)*s*s);
        for (uint i = 0; i < s; ++i) {
                for (uint j = 0; j < s; ++j) {
                        for (uint k = 0; k < s; ++k) {
                                Hnew[j*s + i] += 1e6*Hcm[k*s + j]*Hcm[i*s + k];
                        }
                }
        }
        memcpy(Hcm, Hnew, sizeof(*Hnew)*s*s);
        util_free(Hnew);
}

void example_center_of_mass(double gdd) {
        bool even = true;
        FockBasis basis = fock_basis_alloc_special(2, 41.0, fock_hermite_parity, &even);
        Workspace w = workspace_create_from_basis(&basis);

        // FIXME: This cannot be fixed, maybe add filenames as arguments
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vc_50.dat");
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vdd_50_eta10.dat");

        const uint s = w.basis.size;
        double* Hcm = util_malloc(sizeof(*Hcm)*s*s);
        memset(Hcm, 0, sizeof(*Hcm)*s*s);
        for (uint i = 0; i < s; ++i) {
                Hcm[i*s + i] = fock_energy(&w.basis.states[i]) / 2;
                for (uint j = 0; j < s; ++j) {
                        double v = __Hcm_mval(&w.basis.states[j], &w.basis.states[i]);
                        Hcm[i*s + j] += v;
                }
        }

        __Hcm_square(Hcm, s);
        for (uint i = 0; i < s; ++i)
                Hcm[i*s + i] -= 0.5;

        char fname[256];
        sprintf(fname, "relative_energies_gdd%.1lf.txt", gdd);
        FILE* f = util_fopen(fname, "w");

        double cc[2] = {0.0, gdd};
        for (double g = -50.0; g <= 50.1; g += 0.1) {
                cc[0] = g;
                DEigenSol sol = workspace_compute_Hcm(&w, cc, Hcm);

                fprintf(f, "%lf", g);
                for (uint i = 0; i < sol.size; ++i)
                        fprintf(f, " %e", sol.val[i]);
                putc('\n', f);

                deigen_sol_free(&sol);
        }

        util_fclose(f);
}

void example_center_of_mass_over_a(double gdd) {
        bool even = true;
        FockBasis basis = fock_basis_alloc_special(2, 40.0, fock_hermite_parity, &even);
        Workspace w = workspace_create_from_basis(&basis);

        // FIXME: This cannot be fixed, maybe add filenames as arguments
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vc_50.dat");
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vdd_50_eta10.dat");

        const uint s = w.basis.size;
        double* Hcm = util_malloc(sizeof(*Hcm)*s*s);
        memset(Hcm, 0, sizeof(*Hcm)*s*s);
        for (uint i = 0; i < s; ++i) {
                Hcm[i*s + i] = fock_energy(&w.basis.states[i]) / 2;
                for (uint j = 0; j < s; ++j) {
                        double v = __Hcm_mval(&w.basis.states[j], &w.basis.states[i]);
                        Hcm[i*s + j] += v;
                }
        }

        __Hcm_square(Hcm, s);
        for (uint i = 0; i < s; ++i)
                Hcm[i*s + i] -= 0.5;

        char fname[256];
        sprintf(fname, "./data/energies/relative_energies_over_a_gdd%.1lf.txt", gdd);
        FILE* f = util_fopen(fname, "w");
        fprintf(f, "#   1/a   g    E[0]   E[1]  ...\n");

        const double eta = 10.0;
        const double C = 1.46035450880958681;

        double cc[2] = {0.0, gdd};
        for (double oa = -100.0; oa <= 100.0; oa += 0.2) {
                /* First convert to harmonic oscillator units */
                /*double oa = over_a * 0.05918;*/
                if (oa == C*sqrt(eta))
                        continue;
                double g = 2.0*eta/(oa - C*sqrt(eta));
                cc[0] = g;
                DEigenSol sol = workspace_compute_Hcm(&w, cc, Hcm);

                fprintf(f, "%lf %lf", oa, g);
                // Sentinel for easier data analysis
                int diff = 0;
                if (oa > C*sqrt(eta)) {
                        fprintf(f, " %lf", -2e5);
                        diff = -1;
                }
                for (uint i = 0; i < sol.size+diff; ++i)
                        fprintf(f, " %e", sol.val[i]);
                putc('\n', f);

                deigen_sol_free(&sol);
        }

        util_fclose(f);
}

void example_convergence_over_a(double oa) {
        bool even = true;

        char fname[256];
        sprintf(fname, "./data/convergence_oa%.1lf.txt", oa);
        FILE* f_convergence = util_fopen(fname, "w");

        for (double cutoff = 10.0; cutoff <= 50.0; cutoff += 1.0) {

                FockBasis basis = fock_basis_alloc_special(2, cutoff, fock_hermite_parity, &even);
                Workspace w = workspace_create_from_basis(&basis);

                workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vc_50.dat");

                const uint s = w.basis.size;
                double* Hcm = util_malloc(sizeof(*Hcm)*s*s);
                memset(Hcm, 0, sizeof(*Hcm)*s*s);
                for (uint i = 0; i < s; ++i) {
                        Hcm[i*s + i] = fock_energy(&w.basis.states[i]) / 2;
                        for (uint j = 0; j < s; ++j) {
                                double v = __Hcm_mval(&w.basis.states[j], &w.basis.states[i]);
                                Hcm[i*s + j] += v;
                        }
                }

                __Hcm_square(Hcm, s);
                for (uint i = 0; i < s; ++i)
                        Hcm[i*s + i] -= 0.5;


                const double eta = 10.0;
                const double C = 1.46035450880958681;

                if (oa == C*sqrt(eta))
                        continue;
                double g = 2.0*eta/(oa - C*sqrt(eta));
                DEigenSol sol = workspace_compute_Hcm(&w, &g, Hcm);

                fprintf(f_convergence, "%lf ", cutoff);
                // Sentinel for easier data analysis
                int diff = 0;
                if (oa > C*sqrt(eta)) {
                        fprintf(f_convergence, " %lf", -2e5);
                        diff = -1;
                }
                for (uint i = 0; i < 6+diff; ++i)
                        fprintf(f_convergence, " %e", sol.val[i]);
                putc('\n', f_convergence);

                deigen_sol_free(&sol);
                workspace_delete(&w);
        }

        util_fclose(f_convergence);
}

void example_center_of_mass_print_wf(double g, double gdd) {
        bool even = true;
        FockBasis basis = fock_basis_alloc_special(2, 41.0, fock_hermite_parity, &even);
        Workspace w = workspace_create_from_basis(&basis);

        // FIXME: This cannot be fixed, maybe add filenames as arguments
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vc_50.dat");
        workspace_add_potential(&w, "/home/mateusz/projects/dipoles/data/Vdd_50_eta10.dat");

        const uint s = w.basis.size;
        double* Hcm = util_malloc(sizeof(*Hcm)*s*s);
        memset(Hcm, 0, sizeof(*Hcm)*s*s);
        for (uint i = 0; i < s; ++i) {
                Hcm[i*s + i] = fock_energy(&w.basis.states[i]) / 2;
                for (uint j = 0; j < s; ++j) {
                        double v = __Hcm_mval(&w.basis.states[j], &w.basis.states[i]);
                        Hcm[i*s + j] += v;
                }
        }

        __Hcm_square(Hcm, s);
        for (uint i = 0; i < s; ++i)
                Hcm[i*s + i] -= 0.5;

        char fname[256];
        sprintf(fname, "./data/wf_g%.1lf_gdd%.1lf.txt", g, gdd);
        FILE* f = util_fopen(fname, "w");

        double cc[2] = {g, gdd};
        DEigenSol sol = workspace_compute_Hcm(&w, cc, Hcm);

        for (double x = 0.0; x < 5.0001; x += 0.01) {
                double x12[2] = { 0.5 * x, -0.5 * x};
                fprintf(f, "%lf %e\n", x, fock_basis_compute(&w.basis, sol.U, x12));
        }

        deigen_sol_free(&sol);
        util_fclose(f);
        workspace_delete(&w);
}
