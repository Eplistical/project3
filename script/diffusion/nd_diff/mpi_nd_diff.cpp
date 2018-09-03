#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <iterator>
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fmtstring.hpp"
#include "misc/sparsemat.hpp"
#include "misc/matrixop_io.hpp"
#include "misc/label_index_convertor.hpp"
#include "misc/MPIer.hpp"
#include "ode.hpp"

using namespace std;

using BOOL = bool;
using INTEGER = long int;
using REAL = double;

const INTEGER DIM(3);
const vector<INTEGER> Nx(DIM, 15);
const vector<REAL> dx(DIM, 0.02);
const vector<REAL> D(DIM, 1.0);
const REAL dt(0.00001);
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

void cal_dudt(
        const int* /* NEQ */, 
        const REAL* /* t */, 
        const REAL* u, 
        REAL* dudt,
        double* /* rpar */,
        int* /* ipar */
        )
{
    INTEGER mid, idx_l, idx_r, idx_wall;

    for (INTEGER idx(0); idx < my_Ntot; ++idx) {
        mid = DIM + (idx + my_idx_ofs) * maxNneigh;
        dudt[idx] = 0.0;

        for (INTEGER d(0); d < DIM; ++d) {
            idx_l = neigh_list[mid - d - 1];
            if (idx_l != -1) idx_l -= my_idx_ofs;
            idx_r = neigh_list[mid + d + 1];
            if (idx_r != -1) idx_r -= my_idx_ofs;

            idx_wall = idx % N_per_layer;

            if ( (d == DIM - 1) and ( (idx < N_per_layer) or ( idx > (my_Ntot - N_per_layer) ) ) ) {
                if (idx < N_per_layer) {
                    // first layer, last dimension
                    if (not u_l.empty()) {
                        dudt[idx] += (u_l[idx_wall] + u[idx_r] - 2 * u[idx]) * Dinvdx2[d];
                    }
                    else {
                        dudt[idx] += (u[idx_r] - u[idx]) * Dinvdx2[d];
                    }
                }
                else {
                    // last layer, last dimension
                    if (not u_r.empty()) {
                        dudt[idx] += (u_r[idx_wall] + u[idx_l] - 2 * u[idx]) * Dinvdx2[d];
                    }
                    else {
                        dudt[idx] += (u[idx_l] - u[idx]) * Dinvdx2[d];
                    }
                }
            }
            else {
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
    }
}

void update_ulur(const vector<REAL>& u) {
    // update u_r & u_l values
    MPIer::Request rs_l, rs_r, rr_l, rr_r;
    vector<REAL> send_layer_l, send_layer_r;
    if (MPIer::size > 1) {
        // assign info to send
        if (MPIer::rank == 0) {
            send_layer_l.clear();
            send_layer_r.assign(u.end() - N_per_layer, u.end());
        }
        else if (MPIer::rank == MPIer::size - 1) {
            send_layer_l.assign(u.begin(), u.begin() + N_per_layer);
            send_layer_r.clear();
        }
        else {
            send_layer_l.assign(u.begin(), u.begin() + N_per_layer);
            send_layer_r.assign(u.end() - N_per_layer, u.end());
        }

        // send/recv
        if (MPIer::rank == 0) {
            MPIer::isend(MPIer::rank + 1, send_layer_r, rs_r);
            MPIer::irecv(MPIer::rank + 1, u_r, rr_r);
            MPIer::wait(rs_r, rr_r);
        }
        else if (MPIer::rank == MPIer::size - 1) {
            MPIer::isend(MPIer::rank - 1, send_layer_l, rs_l);
            MPIer::irecv(MPIer::rank - 1, u_l, rr_l);
            MPIer::wait(rs_l, rr_l);
        }
        else {
            MPIer::isend(MPIer::rank - 1, send_layer_l, rs_l);
            MPIer::isend(MPIer::rank + 1, send_layer_r, rs_r);
            MPIer::irecv(MPIer::rank - 1, u_l, rr_l);
            MPIer::irecv(MPIer::rank + 1, u_r, rr_r);
            MPIer::wait(rs_l, rr_l, rs_r, rr_r);
        }
    }
}

int main(int argc, char** argv) {
    MPIer::setup();
    misc::crasher::confirm<>(MPIer::size > 1, "MPI size must > 1!");

    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx

    // work batch
    my_batch = MPIer::assign_job(Nx[DIM-1]);
    my_Ntot = my_batch.size() * N_per_layer;
    my_idx_ofs = my_batch[0] * N_per_layer;

    // init
    init_ijk_list();
    init_neigh_list();
    vector<REAL> u = init_u();

    const REAL atol(1e-8), rtol(1e-3);

    // output header
    if (MPIer::master) {
        out.info("# Nx = ", Nx);
        out.info(misc::fmtstring("# rtol = %.2e, atol = %.2e", rtol, atol));
    }

    const INTEGER Nstep(atoi(argv[1]));
    const INTEGER Anastep(atoi(argv[2]));
    if (MPIer::master) {
        misc::crasher::confirm<>(argc >= 3, "insufficient input para!");
    }

    vector<REAL> u_buf, u_all;

    // loop
    REAL t(0.0), tout(0.0);
    for (INTEGER istep(0); istep < Nstep; ++istep) {
        // output
        if (istep % Anastep == 0) {
            if (MPIer::master) {
                u_all.assign(u.begin(), u.end());
            }

            for (int r(1); r < MPIer::size; ++r) {
                if (MPIer::rank == r) {
                    MPIer::send(0, u);
                }
                else if (MPIer::master) {
                    MPIer::recv(r, u_buf);
                    u_all.insert(u_all.end(), u_buf.begin(), u_buf.end());
                }
                MPIer::barrier();
            }

            if (MPIer::master) {
                out.open(misc::fmtstring("mpind_%d.dat", istep / Anastep), ios::out);
                show_u(u_all, out);
                out.close();
            }
        }

        if (MPIer::master) {
            timer::tic();
        }

        // integrate t -> t+dt
        update_ulur(u);
        tout = t + dt;
        //dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);
        //rkf45(u, t, tout, cal_dudt, rtol, atol);
        //euler(u, t, tout, cal_dudt);
        dopri5(u, t, tout, cal_dudt, rtol, atol);

        // boundary condition
        for (INTEGER idx(0); idx < my_Ntot; ++idx) {
            if (on_the_boundary(idx + my_idx_ofs)) {
                u[idx] = 0.0;
            }
        }

        if (MPIer::master) {
            ioer::tabout("# step ", istep, " done. ", timer::toc());
        }
    }

    MPIer::finalize();
    return 0;

}
