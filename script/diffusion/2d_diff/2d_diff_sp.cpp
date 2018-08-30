#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/fmtstring.hpp"
#include "dvode_f90/dvode.hpp"

using namespace std;

const int Nx(40);
const double dx(0.05);
const double D(1.0);
const double dt(0.0001);

const double Dinvdx2(D / dx / dx);
const int Ntot(Nx*Nx);


vector<double> init_u()
{
    vector<double> u(Ntot, 0.0);
    double mu = 0.5*dx*Nx;
    double sigma = 0.05;

    double x, y;
    int idx;
    double C(pow(2*M_PI * sigma * sigma, -1));
    for (int ix(1); ix < Nx - 1; ++ix) {
        x = ix * dx;
        for (int iy(1); iy < Nx - 1; ++iy) {
            y = iy * dx;
            idx = ix + Nx * iy;
            u[idx] = C * exp(-0.5 * (pow((x-mu) / sigma, 2) + pow((y-mu) / sigma, 2)));
        }
    }
    return u;
}

void cal_dudt(const int* /* NEQ */, const double* /* t */, const double* u, double* dudt)
{
    // calculate dudt
    int idx;

    dudt[0+0*Nx]  = Dinvdx2 * (u[0+1*Nx] + u[1+0*Nx] - 2*u[0+0*Nx]);
    dudt[Nx-1+0*Nx] = Dinvdx2 * (u[Nx-2+0*Nx] + u[Nx-1+1*Nx] - 2*u[Nx-1+0*Nx]);
    dudt[0+(Nx-1)*Nx] = Dinvdx2 * (u[1+(Nx-1)*Nx] + u[0+(Nx-2)*Nx] - 2*u[0+(Nx-1)*Nx]);
    dudt[Nx-1+(Nx-1)*Nx] = Dinvdx2 * (u[Nx-1+(Nx-2)*Nx] + u[(Nx-2)+(Nx-1)*Nx] - 2*u[(Nx-1)+(Nx-1)*Nx]);

    for (int i(1); i < Nx - 1; ++i) {
        dudt[0+i*Nx] = Dinvdx2 * (u[0+(i-1)*Nx] + u[0+(i+1)*Nx] + u[1+i*Nx] - 3*u[0+i*Nx]);
        dudt[Nx-1+i*Nx] = Dinvdx2 * (u[Nx-1+(i-1)*Nx] + u[Nx-1+(i+1)*Nx] + u[Nx-2+i*Nx] - 3*u[Nx-1+i*Nx]);
        dudt[i+0*Nx] = Dinvdx2 * (u[i-1+0*Nx] + u[i+1+0*Nx] + u[i+1*Nx] - 3*u[i+0*Nx]);
        dudt[i+(Nx-1)*Nx] = Dinvdx2 * (u[i-1+(Nx-1)*Nx] + u[i+1+(Nx-1)*Nx] + u[i+(Nx-2)*Nx] - 3*u[i+0*Nx]);
    }

    for (int ix(1); ix < Nx - 1; ++ix) {
        for (int iy(1); iy < Nx - 1; ++iy) {
            idx = ix + Nx * iy;
            dudt[idx] = Dinvdx2 * (u[idx+1] + u[idx-1] + u[idx+Nx] + u[idx-Nx] - 4*u[idx]); 
        }
    }
}

void cal_jac(const int* /* NEQ */, const double* /* t */, const double* u, int* IA, int* JA, int* NZ, double* A)
{
    // calculate j^th column of jacobian matrix
    if (*NZ == 0) {
        *NZ = 5*Ntot*Ntot + 6*Ntot -15;
    }
    else {
        const int ofs(1);
        int count(0);
        IA[0] = ofs;
        int idx;
        for (int ix(0); ix < Nx; ++ix) {
            for (int iy(0); iy < Nx; ++iy) {
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
    ioer::output_t out("");
    out.set_precision(10);

    // init
    int Nstep(atoi(argv[1]));
    double t, tout(dt);
    vector<double> u = init_u();
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");

    // loop
    for (int istep(0); istep < Nstep; ++istep) {
        // output
        out.tabout("# step = ", istep);
        tout = t + dt;

        dvode_sp(u, t, tout, cal_dudt, cal_jac, 1e-3, 1e-14);

        // apply boundary condition
        for (int i(0); i < Nx; ++i) {
            u[i + 0*Nx] = 0.0;
            u[i + (Nx-1)*Nx] = 0.0;
        }
        for (int i(0); i < Nx; ++i) {
            u[0 + i*Nx] = 0.0;
            u[Nx-1 + i*Nx] = 0.0;
        }
    }

    for (int ix(0); ix < Nx; ++ix) {
        for (int iy(0); iy < Nx; ++iy) {
            printf("%16.6f%16.6f%16.6f\n", ix*dx, iy*dx, u[ix + iy*Nx]);
        }
    }

    return 0;
}
