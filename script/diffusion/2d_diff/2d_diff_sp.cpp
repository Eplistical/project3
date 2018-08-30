#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fmtstring.hpp"
#include "misc/matrixop_io.hpp"
#include "dvode.hpp"

using namespace std;

const int Nx(40);
const double dx(0.05);
const double D(1.0);
const double dt(0.0001);

const double Dinvdx2(D / dx / dx);
const int Ntot(Nx*Nx);

inline int get_idx(const int ix, const int iy) {
    // idx -> ijk
    if (ix < 0 or ix >= Nx) return -1;
    if (iy < 0 or iy >= Nx) return -1;
    return ix + iy*Nx;
}

inline void get_ijk(const int idx, int& ix, int& iy) {
    // ijk -> idx
    ix = idx % Nx;
    iy = idx / Nx;
}

inline bool on_the_boundary(const int idx) {
    // true if the point is on the boundary
    int ix, iy;
    get_ijk(idx, ix, iy);
    if (ix == 0 or ix == Nx-1) return true;
    if (iy == 0 or iy == Nx-1) return true;
    return false;
}

inline void show_u(vector<double>& u, ioer::output_t& out = ioer::STDOUT) {
    // output u
    int ix, iy;
    for (int idx(0); idx < Ntot; ++idx) {
        get_ijk(idx, ix, iy);
        out.tabout(ix*dx, iy*dx, misc::fmtstring(" %20.6e", u[idx]));
    }
}

inline void find_neigh(const int idx, vector<int>& neigh_idx) {
    // find neighbor indices (accending order) for a given point index, 
    //  the list includes the point itself

    int ix, iy;
    get_ijk(idx, ix, iy);

    int N_neigh(1);
    N_neigh += (ix == 0 or ix == Nx - 1) ? 1 : 2;
    N_neigh += (iy == 0 or iy == Nx - 1) ? 1 : 2;

    neigh_idx.clear();
    neigh_idx.reserve(N_neigh);

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
}

vector<double> init_u()
{
    // initialize u
    vector<double> u(Ntot, 0.0);
    double mu = 0.5*dx*(Nx);
    double sigma = 0.05;

    double x, y;
    double C(pow(2*M_PI * sigma * sigma, -1));
    for (int ix(0); ix < Nx; ++ix) {
        x = ix * dx;
        for (int iy(0); iy < Nx; ++iy) {
            y = iy * dx;

            u[get_idx(ix,iy)] = C * exp(-0.5 * (pow((x-mu) / sigma, 2) + pow((y-mu) / sigma, 2)));
        }
    }
    return u;
}

void cal_dudt(const int* /* neq */, const double* /* t */, const double* u, double* dudt)
{
    // calculate dudt
    int idx, idx2;
    int Nneigh(0);

    for (int ix(0); ix < Nx; ++ix) {
        for (int iy(0); iy < Nx; ++iy) {
            idx = get_idx(ix, iy);
            dudt[idx] = 0.0;

            Nneigh = 0;
            idx2 = get_idx(ix-1, iy);
            if (idx2 != -1) {
                dudt[idx] += u[idx2];
                Nneigh += 1;
            }

            idx2 = get_idx(ix+1, iy);
            if (idx2 != -1) {
                dudt[idx] += u[idx2];
                Nneigh += 1;
            }

            idx2 = get_idx(ix, iy-1);
            if (idx2 != -1) {
                dudt[idx] += u[idx2];
                Nneigh += 1;
            }

            idx2 = get_idx(ix, iy+1);
            if (idx2 != -1) {
                dudt[idx] += u[idx2];
                Nneigh += 1;
            }
            dudt[idx] = Dinvdx2 * (dudt[idx] - Nneigh * u[idx]);
        } 
    }
}

void cal_jac(const int* /* NEQ */, const double* /* t */, const double* u, int* IA, int* JA, int* NZ, double* A)
{
    // calculate j^th column of jacobian matrix
    if (*NZ == 0) {
        *NZ = 5*Nx*Nx + 6*Nx - 15;
    }
    else {
        const int ofs(1);
        int count(0);
        IA[0] = ofs;
        int idx;
        for (int iy(0); iy < Nx; ++iy) {
            for (int ix(0); ix < Nx; ++ix) {
                idx = ix + iy * Nx;

                if (ix == 0) {
                    if (iy == 0) {
                        JA[count] = idx    + ofs;       A[count] = -2*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                    else if (iy == Nx-1) {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -2*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                    else {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -3*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                }
                else if (ix == Nx-1) {
                    if (iy == 0) {
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -2*Dinvdx2;         count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                    else if (iy == Nx-1) {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -2*Dinvdx2;         count += 1;
                    }
                    else {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -3*Dinvdx2;         count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                }
                else {
                    if (iy == 0) {
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -3*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                    else if (iy == Nx-1) {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -3*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                    else {
                        JA[count] = idx-Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx-1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx    + ofs;       A[count] = -4*Dinvdx2;         count += 1;
                        JA[count] = idx+1  + ofs;       A[count] = Dinvdx2;            count += 1;
                        JA[count] = idx+Nx + ofs;       A[count] = Dinvdx2;            count += 1;
                    }
                }
                IA[idx+1] = count + ofs;
            }
        }
    }
}

int main(int argc, char** argv) {
    // output 
    ioer::output_t out;
    out.set_precision(10);

    // mem approx
    ioer::info("# Nx = ", Nx);
    int nz1 = 5*Nx*Nx + 6*Nx - 15;
    double mem = (nz1 + Nx + 1) * 4 + nz1 * 8;
    mem = mem / 1024 / 1024 / 1024;
    ioer::info(misc::fmtstring("# approx mem = %.6f GB", mem));


    // init
    int Nstep(atoi(argv[1]));
    double t, tout(dt);
    vector<double> u = init_u();
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");

    const double atol(1e-12), rtol(1e-4);
    ioer::info(misc::fmtstring("# rtol = %.2e, atol = %.2e", rtol, atol));

    int ix, iy;
    int Anastep(20);

    // loop
    for (int istep(0); istep < Nstep; ++istep) {
        // output
        if (istep % Anastep == 0) {
            out.open(misc::fmtstring("2d_%d.dat", istep / Anastep), ios::out);
            show_u(u, out);
            out.close();
        }

        timer::tic();

        tout = t + dt;

        dvode_sp(u, t, tout, cal_dudt, nullptr, rtol, atol);
        //dvode_sp(u, t, tout, cal_dudt, cal_jac, rtol, atol);

        // apply boundary condition
        for (int idx(0); idx < Ntot; ++idx) {
            if (on_the_boundary(idx)) {
                u[idx] = 0.0;
            }
        }
        out.tabout("# ", timer::toc());
    }

    return 0;
}
