#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fmtstring.hpp"
#include "misc/sparsemat.hpp"
#include "misc/matrixop_io.hpp"
#include "misc/label_index_convertor.hpp"
#include "misc/MPIer.hpp"
#include "dvode.hpp"

using namespace std;

using BOOL = bool;
using INTEGER = long int;
using REAL = double;


const INTEGER DIM(2);
const vector<INTEGER> Nx(DIM, 4);
const vector<REAL> dx(DIM, 0.05);
const vector<REAL> D(DIM, 1.0);
const REAL dt(0.0001);
const vector<REAL> Dinvdx2(D / dx / dx);
const INTEGER Ntot(product(Nx));


const INTEGER Nlayer(Nx[DIM-1]); // # layers
const INTEGER N_per_layer(Ntot / Nlayer); // N points in a single layer 

vector<INTEGER> my_batch;
INTEGER my_Ntot;
INTEGER my_idx_ofs;
vector<REAL> u_l; 
vector<REAL> u_r; 

const INTEGER maxNneigh(2*DIM+1);

vector<INTEGER> neigh_list;
vector<INTEGER> ijk_list; 

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

inline INTEGER get_Nneigh(const INTEGER idx) {
    // get total number of neighbors for given point idx, the result excludes the point itself
    const INTEGER* neigh(&neigh_list[idx * maxNneigh]);
    INTEGER Nneigh(0);
    for (INTEGER k(0); k < maxNneigh; ++k) {
        if (neigh[k] != -1 and k != DIM)
            Nneigh += 1;
    }
    return Nneigh;
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
    const vector<REAL> mu(0.5 * dx * Nx);
    const vector<REAL> sigma(DIM, 0.05);
    const REAL C(1.0 / pow(2 * M_PI, 0.5 * DIM) / product(sigma));

    vector<REAL> u(my_Ntot, 0.0);

    vector<INTEGER> ijk;
    vector<REAL> xyz;
    for (INTEGER idx(0); idx < my_Ntot; ++idx) {
        get_ijk(my_idx_ofs + idx, ijk);
        xyz = ijk * dx;

        u[idx] = C * exp(-0.5 * sum(pow((xyz - mu) / sigma, 2)));
    }
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

#ifdef ENABLE_TEST
    ioer::newline();
    ioer::info(" # -- TEST FOR cal_dudt BEGIN");
    ioer::tabout("idx", "u[idx]", "dudt[idx]");
    for (INTEGER idx(0); idx < Ntot; ++idx) {
        ioer::tabout(idx, u[idx], dudt[idx]);
    }
    ioer::info(" # -- TEST FOR cal_dudt END");
    ioer::newline();
#endif

}


void cal_jac(const int* /* NEQ */, const REAL* /* t */, const REAL* u, int* IA, int* JA, int* NZ, REAL* A)
{
    // calculate jacobian matrix in sparse form
    if (*NZ == 0) {
        // if *NZ is 0, return number of non-zero elements + diag elements
        INTEGER NZ_tmp(0);
        for (INTEGER idx(0); idx < Ntot; ++idx) {
            NZ_tmp += get_Nneigh(idx);
        }
        NZ_tmp += Ntot;

        misc::crasher::confirm<overflow_error>(NZ_tmp <= numeric_limits<int>::max(), "NZ_tmp > max int, overflow");
        *NZ = static_cast<int>(NZ_tmp);
    }
    else {
        // if *NZ is not 0, calculate sparse Jacobian matrix in Fortran style
        const int ofs(1);
        INTEGER count(0);
        INTEGER mid;
        INTEGER Nneigh;
        REAL A_self;
        INTEGER A_self_count;
        IA[0] = ofs;

        for (INTEGER idx(0); idx < Ntot; ++idx) {
            Nneigh = get_Nneigh(idx);
            mid = idx * maxNneigh + DIM;
            A_self = 0.0;

            for (INTEGER k(-DIM); k <= DIM; ++k) {
                if (neigh_list[mid + k] != -1) {
                    JA[count] = neigh_list[mid + k] + ofs;
                    if (k == 0) {
                        A_self_count = count;
                    }
                    else {
                        A[count] = Dinvdx2[abs(k) - 1];
                        A_self -= Dinvdx2[abs(k) - 1];
                    }
                    count += 1;
                }
            }
            A[A_self_count] = A_self;
            IA[idx + 1] = count + ofs;
        }

#ifdef ENABLE_TEST
    ioer::newline();
    ioer::info(" # -- TEST FOR cal_jac BEGIN");
    matrixop::showvec("u", u, Ntot);
    matrixop::showvec("A", A, count);
    matrixop::showvec("JA", JA, count);
    matrixop::showvec("IA", IA, Ntot);
    ioer::info(" # -- TEST FOR cal_jac END");
    ioer::newline();
#endif

    }
}

int main(int argc, char** argv) {
    MPIer::setup();
    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx
    out.info("# Nx = ", Nx);

    // work batch
    my_batch = MPIer::assign_job(Nx[DIM-1]);
    my_Ntot = my_batch.size() * N_per_layer;
    my_idx_ofs = my_batch[0] * N_per_layer;

    // init
    init_ijk_list();
    init_neigh_list();
    vector<REAL> u = init_u();
    u_l.resize(N_per_layer);
    u_r.resize(N_per_layer);

    //misc::crasher::confirm<>(argc >= 2, "insufficient input para!");

    const REAL atol(1e-8), rtol(1e-3);
    ioer::info(misc::fmtstring("# rtol = %.2e, atol = %.2e", rtol, atol));

    //INTEGER Nstep(atoi(argv[1]));
    const INTEGER Nstep(21);
    const INTEGER Anastep(20);

    // loop
    REAL t(0.0), tout(0.0);
    for (INTEGER istep(0); istep < Nstep; ++istep) {
        timer::tic();

        // update u_r & u_l values
        MPIer::Request rs_l, rs_r, rr_l, rr_r;
        vector<double> send_layer_l, send_layer_r;
        if (MPIer::size > 1) {
            send_layer_l.assign(u.begin(), u.begin() + N_per_layer);
            send_layer_r.assign(u.end() - N_per_layer, u.end());

            if (MPIer::rank == 0) {
                MPIer::isend(MPIer::rank + 1, send_layer_r, rs_r);
                MPIer::irecv(MPIer::rank + 1, u_r, rr_r);
                rs_r.Wait();
                rr_r.Wait();
            }
            else if (MPIer::rank == MPIer::size - 1) {
                MPIer::isend(MPIer::rank - 1, send_layer_l, rs_l);
                MPIer::irecv(MPIer::rank - 1, u_l, rr_l);
                rs_l.Wait();
                rr_l.Wait();
            }
            else {
                MPIer::isend(MPIer::rank - 1, send_layer_l, rs_l);
                MPIer::isend(MPIer::rank + 1, send_layer_r, rs_r);
                MPIer::irecv(MPIer::rank - 1, u_l, rr_l);
                MPIer::irecv(MPIer::rank + 1, u_r, rr_r);
                rs_l.Wait();
                rs_r.Wait();
                rr_l.Wait();
                rr_r.Wait();
            }
        }
        else {
        }
        ioer::info("rank ", MPIer::rank, " my_batch = ", my_batch);
        ioer::info("rank ", MPIer::rank, " my_Ntot = ", my_Ntot);
        ioer::info("rank ", MPIer::rank, " my_idx_ofs = ", my_idx_ofs);
        ioer::info("rank ", MPIer::rank, " u = ", u);
        ioer::info("rank ", MPIer::rank, " u_l = ", u_l);
        ioer::info("rank ", MPIer::rank, " u_r = ", u_r);
        break;

        // integrate
        tout = t + dt;
        dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);

        // boundary condition
        for (INTEGER idx(0); idx < my_Ntot; ++idx) {
            if (on_the_boundary(idx + my_idx_ofs)) {
                u[idx] = 0.0;
            }
        }
        ioer::tabout("# step ", istep, " done. ", timer::toc());
    }

    MPIer::finalize();
    return 0;

    /*
    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx
    out.info("# Nx = ", Nx);
    INTEGER nz1 = 7*Nx*Nx*Nx - 6*Nx*Nx;
    REAL mem = (nz1 + Nx + 1) * 4 + nz1 * 8 + (Ntot * 7) * 4;
    mem = mem / 1024 / 1024 / 1024;
    out.info(misc::fmtstring("# approx mem = %.6f GB", mem));

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
    */
    for (int istep(0); istep < Nstep; ++istep) {
        if (istep % Anastep == 0) {
            out.open(misc::fmtstring("%d.dat", istep / Anastep), ios::out);
            show_u(u, out);
            out.close();
        }

        // output
        timer::tic();

        tout = t + dt;

        //dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);
        dvode_sp(u, t, tout, cal_dudt, cal_jac, rtol, atol);

        // apply boundary condition
        for (INTEGER idx(0); idx < Ntot; ++idx) {
            if (on_the_boundary(idx)) {
                u[idx] = 0.0;
            }
        }
        ioer::tabout("# step ", istep, " done. ", timer::toc());
    }

    MPIer::finalize();
    return 0;
}
