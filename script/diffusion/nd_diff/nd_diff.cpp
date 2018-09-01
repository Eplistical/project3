#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fmtstring.hpp"
#include "misc/sparsemat.hpp"
#include "misc/matrixop_io.hpp"
#include "misc/label_index_convertor.hpp"
#include "dvode.hpp"

using namespace std;

using BOOL = bool;
using INTEGER = long int;
using REAL = double;


const INTEGER DIM(3);
const vector<INTEGER> Nx(DIM, 15);
const vector<REAL> dx(DIM, 0.02);
const vector<REAL> D(DIM, 1.0);
const REAL dt(0.0001);
const INTEGER Ntot(product(Nx));
const vector<REAL> Dinvdx2(D / dx / dx);

const INTEGER maxNneigh(2*DIM+1);

vector<INTEGER> neigh_list;
vector<INTEGER>ijk_list; 

inline void get_ijk(const INTEGER idx, vector<INTEGER>& ijk) {
    ijk.assign(ijk_list.begin() + idx * DIM, ijk_list.begin() + (idx + 1) * DIM);
}

inline bool on_the_boundary(const INTEGER idx) {
    INTEGER mid(DIM + idx * maxNneigh);
    for (INTEGER d(0); d < DIM; ++d) {
        if (neigh_list[mid + d + 1] == -1) return true;
        if (neigh_list[mid - d - 1] == -1) return true;
    }
    return false;
}

void show_u(vector<REAL>& u, ioer::output_t& out = ioer::STDOUT) {
    // output u
    vector<INTEGER> ijk;
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ijk);
        out.tabout(ijk, ijk*dx, misc::fmtstring(" %20.6e", u[idx]));
    }
}

void init_ijk_list() {
    misc::LabelIndexConverter<INTEGER> cnv(Nx);
    vector<INTEGER> ijk(DIM, 0);
    ijk_list.reserve(Ntot * DIM);

    for (INTEGER idx(0); idx < Ntot; ++idx) {
        ijk = cnv.LabelIndexConverter::index_to_label(idx);
        ijk_list.insert(ijk_list.end(), ijk.begin(), ijk.end());
    }

#ifdef ENABLE_TEST
    ioer::newline();
    ioer::info(" # -- TEST FOR init_ijk_list BEGIN");
    ioer::tabout("idx", "ijk");
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        ioer::tabout(idx, vector<INTEGER>(ijk_list.begin() + DIM * idx, ijk_list.begin() + DIM * (idx + 1)));
    }
    ioer::info(" # -- TEST FOR init_ijk_list END");
    ioer::newline();
#endif

}

void init_neigh_list() {
    neigh_list.resize(Ntot * maxNneigh, -1);
    vector<INTEGER> ijk;

    for (INTEGER idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ijk);

        INTEGER mid(DIM + idx * maxNneigh);

        neigh_list[mid] = idx;

        INTEGER factor(1);
        for (INTEGER d(0); d < DIM; ++d) {
            if (ijk[d] > 0) {
                neigh_list[mid - d - 1] = idx - factor;
            }
            if (ijk[d] < Nx[d]-1) {
                neigh_list[mid + d + 1] = idx + factor;
            }
            factor *= Nx[d];
        }
    }

#ifdef ENABLE_TEST
    ioer::newline();
    ioer::info(" # -- TEST FOR init_neigh_list BEGIN");
    ioer::tabout("idx", "on the boundary", "neigh_list");
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        ioer::tabout(idx, on_the_boundary(idx), vector<INTEGER>(neigh_list.begin() + idx * maxNneigh, neigh_list.begin() + (idx + 1) * maxNneigh));
    }
    ioer::info(" # -- TEST FOR init_negih_list END");
    ioer::newline();
#endif

}

vector<REAL> init_u()
{
    // initialize u
    vector<REAL> u(Ntot, 0.0);
    vector<REAL> mu(0.5*dx*Nx);
    vector<REAL> sigma(DIM, 0.05);

    const REAL C(1.0 / pow(2 * M_PI, 0.5 * DIM) / product(sigma));

    vector<INTEGER> ijk;
    vector<REAL> xyz;
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ijk);
        xyz = ijk * dx;

        u[idx] = C * exp(-0.5 * sum(pow((xyz - mu) / sigma, 2)));
    }

#ifdef ENABLE_TEST
    ioer::newline();
    ioer::info(" # -- TEST FOR init_u BEGIN");
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ijk);
        xyz = ijk * dx;
        ioer::tabout(xyz, u[idx]);
    }
    ioer::info(" # -- TEST FOR init_u END");
    ioer::newline();
#endif

    return u;
}

void cal_dudt(const int* /* NEQ */, const REAL* /* t */, const REAL* u, REAL* dudt)
{
    // calculate dudt
    INTEGER mid, idx_l, idx_r;
    REAL uidx;
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        mid = DIM + idx * maxNneigh;

        dudt[idx] = 0.0;
        for (INTEGER d(0); d < DIM; ++d) {
            idx_l = neigh_list[mid - d - 1];
            idx_r = neigh_list[mid + d + 1];

            if (idx_l == -1) {
               dudt[idx] += (u[idx_r] - u[idx]) * Dinvdx2[d];
            }
            else if (idx_r == -1) {
               dudt[idx] += (u[idx_l] - u[idx]) * Dinvdx2[d];
            }
            else {
               dudt[idx] += (u[idx_l] + u[idx_r] - 2 * u[idx]) * Dinvdx2[d];
            }
        }
    }
}


#if 0
void cal_jac(const int* /* NEQ */, const REAL* /* t */, const REAL* u, int* IA, int* JA, int* NZ, REAL* A)
{
    // calculate jacobian matrix in sparse form
    if (*NZ == 0) {
        *NZ = 7*Nx*Nx*Nx - 6*Nx*Nx;
    }
    else {
        const int ofs(1);
        int count(0);
        int idx;
        int Nneigh;
        IA[0] = ofs;

        for (int idx(0); idx < Ntot; ++idx) {
            vector<int>& neigh = neigh_list.at(idx);
            Nneigh = neigh.size();

            for (int k(0); k < Nneigh; ++k) {
                JA[count] = neigh[k] + ofs;
                if (neigh[k] == idx) {
                    A[count] = -Dinvdx2 * (Nneigh - 1);
                }
                else {
                    A[count] = Dinvdx2;
                }
                count += 1;
            }
            IA[idx+1] = count + ofs;
        }
    }

}
#endif

int main(int argc, char** argv) {
    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx
    out.info("# Nx = ", Nx);
    /*
    INTEGER nz1 = 7*Nx*Nx*Nx - 6*Nx*Nx;
    REAL mem = (nz1 + Nx + 1) * 4 + nz1 * 8 + (Ntot * 7) * 4;
    mem = mem / 1024 / 1024 / 1024;
    out.info(misc::fmtstring("# approx mem = %.6f GB", mem));
    */

    // init
    init_ijk_list();
    init_neigh_list();
    vector<REAL> u = init_u();
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");


    const REAL atol(1e-8), rtol(1e-3);
    ioer::info(misc::fmtstring("# rtol = %.2e, atol = %.2e", rtol, atol));

    INTEGER Nstep(atoi(argv[1]));
    INTEGER Anastep(20);

    // loop
    REAL t(0.0), tout(0.0);
    for (int istep(0); istep < Nstep; ++istep) {
        if (istep % Anastep == 0) {
            out.open(misc::fmtstring("%d.dat", istep / Anastep), ios::out);
            show_u(u, out);
            out.close();
        }

        // output
        timer::tic();

        tout = t + dt;

        dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);
        //dvode_sp(u, t, tout, cal_dudt, cal_jac, rtol, atol);

        // apply boundary condition
        for (INTEGER idx(0); idx < Ntot; ++idx) {
            if (on_the_boundary(idx)) {
                u[idx] = 0.0;
            }
        }
        ioer::tabout("# step ", istep, " done. ", timer::toc());
    }


    return 0;
}
