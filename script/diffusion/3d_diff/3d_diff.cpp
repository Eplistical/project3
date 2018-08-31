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
#include "dvode.hpp"

using namespace std;

const int Nx(15);
const double dx(0.02);
const double D(1.0);
const double dt(0.0001);

const double Dinvdx2(D / dx / dx);

const int Nx2(Nx*Nx);
const int Ntot(Nx*Nx*Nx);

inline int get_idx(const int ix, const int iy, const int iz) {
    // ijk -> idx
    if (ix < 0 or ix >= Nx) return -1;
    if (iy < 0 or iy >= Nx) return -1;
    if (iz < 0 or iz >= Nx) return -1;
    return ix + iy*Nx + iz*Nx2;
}

inline void show_u(vector<double>& u, ioer::output_t& out = ioer::STDOUT) {
    // output u
    for (int iz(0); iz < Nx; ++iz) {
        for (int iy(0); iy < Nx; ++iy) {
            for (int ix(0); ix < Nx; ++ix) {
                out.tabout(ix, iy, iz, ix*dx, iy*dx, iz*dx,
                            misc::fmtstring(" %20.6e", u[get_idx(ix,iy,iz)]));
            }
        }
    }
}

int find_neigh(const int ix, const int iy, const int iz, vector<int>& neigh_idx) 
{
    // find neighbor indices (accending order) for a given point index, 
    //  the list includes the point itself

    int idx = get_idx(ix, iy, iz);
    int N_neigh(1);

    N_neigh += (ix == 0 or ix == Nx - 1) ? 1 : 2;
    N_neigh += (iy == 0 or iy == Nx - 1) ? 1 : 2;
    N_neigh += (iz == 0 or iz == Nx - 1) ? 1 : 2;

    neigh_idx.resize(N_neigh);
    int count(0);

    if (iz != 0) {
        neigh_idx[count] = idx - Nx2;
        count += 1;
    }
    if (iy != 0) {
        neigh_idx[count] = idx - Nx;
        count += 1;
    }
    if (ix != 0) {
        neigh_idx[count] = idx - 1;
        count += 1;
    }

    neigh_idx[count] = idx;
    count += 1;

    if (ix != Nx-1) {
        neigh_idx[count] = idx + 1;
        count += 1;
    }
    if (iy != Nx-1) {
        neigh_idx[count] = idx + Nx;
        count += 1;
    }
    if (iz != Nx-1) {
        neigh_idx[count] = idx + Nx2;
        count += 1;
    }
    return N_neigh;
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
    int Nneigh;

    for (int iz(0); iz < Nx; ++iz) {
        for (int iy(0); iy < Nx; ++iy) {
            for (int ix(0); ix < Nx; ++ix) {
                idx = get_idx(ix, iy, iz);

                Nneigh = 0;
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
        *NZ = 7*Nx*Nx*Nx - 6*Nx*Nx;
    }
    else {
        const int ofs(1);
        int count(0);
        IA[0] = ofs;
        int idx;

        vector<int> neigh;
        int Nneigh;
        for (int iz(0); iz < Nx; ++iz) {
            for (int iy(0); iy < Nx; ++iy) {
                for (int ix(0); ix < Nx; ++ix) {
                    Nneigh = find_neigh(ix, iy, iz, neigh);
                    idx = get_idx(ix, iy, iz);

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
        //dvode_sp(u, t, tout, cal_dudt, cal_jac, rtol, atol);

        // apply boundary condition
        for (int i(0); i < Nx; ++i) {
            for (int j(0); j < Nx; ++j) {
                u[get_idx(0,i,j)]       = 0.0;
                u[get_idx(Nx-1,i,j)]    = 0.0;
                u[get_idx(i,0,j)]       = 0.0;
                u[get_idx(i,Nx-1,j)]    = 0.0;
                u[get_idx(i,j,0)]       = 0.0;
                u[get_idx(i,j,Nx-1)]    = 0.0;
            }
        }
        ioer::tabout("# step ", istep, " done. ", timer::toc());
    }


    return 0;
}
