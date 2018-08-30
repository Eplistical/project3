#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fmtstring.hpp"
#include "dvode.hpp"

using namespace std;

const int Nx(21);
const double dx(0.02);
const double D(0.1);
const double dt(0.000001);

const double Dinvdx2(D / dx / dx);

const int Nx2(Nx*Nx);
const int Ntot(Nx*Nx*Nx);

inline int get_idx(const int ix, const int iy, const int iz) {
    // idx -> ijk
    if (ix < 0 or ix >= Nx) return -1;
    if (iy < 0 or iy >= Nx) return -1;
    if (iz < 0 or iz >= Nx) return -1;
    return ix + iy*Nx + iz*Nx2;
}

inline void get_ijk(const int idx, int& ix, int& iy, int& iz) {
    // ijk -> idx
    ix = idx % Nx;
    iy = (idx % Nx2) / Nx;
    iz = idx / Nx2;
}

inline bool on_the_boundary(const int idx) {
    // true if the point is on the boundary
    int ix, iy, iz;
    get_ijk(idx, ix, iy, iz);
    if (ix == 0 or ix == Nx-1) return true;
    if (iy == 0 or iy == Nx-1) return true;
    if (iz == 0 or iz == Nx-1) return true;
    return false;
}

inline void show_u(vector<double>& u, ioer::output_t& out = ioer::STDOUT) {
    // output u
    int ix, iy, iz;
    for (int idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ix, iy, iz);
        out.tabout(ix, iy, iz, ix*dx, iy*dx, iz*dx, misc::fmtstring(" %20.6e", u[idx]));
    }
}

void find_neigh(const int idx, vector<int>& neigh_idx) {
    // find neighbor indices (accending order) for a given point index, 
    //  the list includes the point itself

    int ix, iy, iz;
    get_ijk(idx, ix, iy, iz);

    int N_neigh(1);
    N_neigh += (ix == 0 or ix == Nx - 1) ? 1 : 2;
    N_neigh += (iy == 0 or iy == Nx - 1) ? 1 : 2;
    N_neigh += (iz == 0 or iz == Nx - 1) ? 1 : 2;

    neigh_idx.clear();
    neigh_idx.reserve(N_neigh);

    if (iz != 0) {
        neigh_idx.push_back(idx - Nx2);
    }
    if (iy != 0) {
        neigh_idx.push_back(idx - Nx);
    }
    if (ix != 0) {
        neigh_idx.push_back(idx - 1);
    }
    neigh_idx.push_back(idx);
    if (ix != Nx-1) {
        neigh_idx.push_back(idx + 1);
    }
    if (iy != Nx-1) {
        neigh_idx.push_back(idx + Nx);
    }
    if (iz != Nx-1) {
        neigh_idx.push_back(idx + Nx2);
    }
}

vector<double> init_u()
{
    // initialize u
    vector<double> u(Ntot, 0.0);
    double mu = 0.5*dx*(Nx-1);
    double sigma = 0.05;

    double x, y, z;
    double C(pow(2*M_PI * sigma * sigma, -1.5));
    for (int ix(0); ix < Nx; ++ix) {
        x = ix * dx;
        for (int iy(0); iy < Nx; ++iy) {
            y = iy * dx;
            for (int iz(0); iz < Nx; ++iz) {
                z = iz * dx;
                u[get_idx(ix,iy,iz)] = C * exp(-0.5 * (pow((x-mu) / sigma, 2) + pow((y-mu) / sigma, 2) + pow((z-mu) / sigma, 2)));
            }
        }
    }
    return u;
}

void cal_dudt(const int* /* NEQ */, const double* /* t */, const double* u, double* dudt)
{
    // calculate dudt
    int idx, idx2;
    int Nneigh(0);

    for (int ix(0); ix < Nx; ++ix) {
        for (int iy(0); iy < Nx; ++iy) {
            for (int iz(0); iz < Nx; ++iz) {
                idx = get_idx(ix, iy, iz);
                dudt[idx] = 0.0;

                idx2 = get_idx(ix-1, iy, iz);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                idx2 = get_idx(ix+1, iy, iz);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                idx2 = get_idx(ix, iy-1, iz);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                idx2 = get_idx(ix, iy+1, iz);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                idx2 = get_idx(ix, iy, iz-1);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                idx2 = get_idx(ix, iy, iz+1);
                if (idx2 != -1) {
                    dudt[idx] += u[idx2] - u[idx];
                    Nneigh += 1;
                }

                dudt[idx] = Dinvdx2 * (dudt[idx] - Nneigh * u[idx]);
            }
        }
    }
}

void cal_jac(const int* /* NEQ */, const double* /* t */, const double* u, int* IA, int* JA, int* NZ, double* A)
{
    // calculate j^th column of jacobian matrix
    if (*NZ == 0) {
     //   *NZ = 7*Nx*Nx*Nx - 6*Nx*Nx;
    }
    else {

    }
}

int main(int argc, char** argv) {
    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx
    out.info("# Nx = ", Nx);
    int nz1 = 7*Nx*Nx*Nx - 6*Nx*Nx;
    double mem = (nz1 + Nx + 1) * 4 + nz1 * 8;
    mem = mem / 1024 / 1024 / 1024;
    out.info(misc::fmtstring("# approx mem = %.6f GB", mem));

    // init
    int Nstep(atoi(argv[1]));
    double t, tout(dt);
    vector<double> u = init_u();
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");

    const double atol(1e-8), rtol(1e-3);
    ioer::info(misc::fmtstring("# rtol = %.2e, atol = %.2e", rtol, atol));

    int ix, iy, iz;
    int Anastep(10);

    // loop
    for (int istep(0); istep < Nstep; ++istep) {
        if (istep % Anastep == 0) {
            out.open(misc::fmtstring("3d_%d.dat", istep / Anastep), ios::out);
            show_u(u, out);
            out.close();
        }

        // output
        timer::tic();

        tout = t + dt;

        dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);
        //dvode_sp(u, t, tout, cal_dudt, cal_jac, 1e-3, 1e-8);

        // apply boundary condition
        for (int idx(0); idx < Ntot; ++idx) {
            if (on_the_boundary(idx)) {
                u[idx] = 0.0;
            }
        }
        ioer::tabout("# step ", istep, " done. ", timer::toc());
    }


    return 0;
}
