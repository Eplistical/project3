#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include <cmath>
#include <cstdlib>
#include <algorithm>

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

extern "C" {
    using cal_dudt_t = void (*)(int* /* NEQ */, 
            double* /* t */, double* /* u */, 
            double* /* dudt */, 
            double* /* rpar */, int* /* ipar */);

    using cal_jac_t = void (*)(int* /* NEQ */, 
            double* /* t */, double* /* u */, 
            int* /* ML */, int* /* MU */, 
            double* /* jac */, int* /* NRPD */, 
            double* /* rpar */, int* /* ipar */);

    extern void dvode_(cal_dudt_t /* cal_dudt */, int* /* NEQ */,  
            double* /* u */, double* /* t */, double* /* tout */,
            int* /* itol */, double* /* rtol */, double* /* atol */,
            int* /* itask */, int* /* istate */, int* /* iopt */, double* /* rwork */,
            int* /* lrw */, int* /* iwork */, int* /* liw */,
            cal_jac_t /* cal_jac */, int* /* MF */, 
            double* /* rpar */, int* /* ipar */
            );
};

void cal_dudt(int* /* NEQ */, double* /* t */, 
        double* u, double* dudt, 
        double* /* rpar */, int* /* ipar */) 
{
    int idx;

    for (int ix(1); ix < Nx - 1; ++ix) {
        for (int iy(1); iy < Nx - 1; ++iy) {
            idx = ix + Nx * iy;
            dudt[idx] = Dinvdx2 * (u[idx+1] + u[idx-1] + u[idx+Nx] + u[idx-Nx] - 4*u[idx]); 
        }
    }
}

void cal_jac(int* /* NEQ */, double* /* t */, 
        double* u, 
        int* /* ML */, int* /* MU */, 
        double* jac,  int* /* NRPD */,
        double* /* rpar */, int* /* ipar */) 
{
    int idx;

    for (int ix(1); ix < Nx - 1; ++ix) {
        for (int iy(1); iy < Nx - 1; ++iy) {
            idx = ix + Nx * iy;
            jac[idx + idx*Ntot] = -4*Dinvdx2;
            jac[idx + (idx+1)*Ntot] = Dinvdx2;
            jac[idx + (idx-1)*Ntot] = Dinvdx2;
            jac[idx + (idx+Nx)*Ntot] = Dinvdx2;
            jac[idx + (idx-Nx)*Ntot] = Dinvdx2;
        }
    }
}

int main(int argc, char** argv) {
    // output 
    ioer::output_t out("");
    out.set_precision(10);

    // integrator para
    int NEQ(Ntot);
    int MF(21); // BDF
    double rtol(1e-4);
    double atol(1e-14);
    int itol(1), itask(1), iopt(0);
    int lrw(22 + 9*NEQ + 2*NEQ*NEQ);
    int liw(30 + NEQ);
    vector<double> rwork(lrw, 0.0);
    vector<int> iwork(liw, 0.0);
    vector<double> rpar(1, 0.0);
    vector<int> ipar(1, 0);
    int istate(1);

    // init
    vector<double> u = init_u();
    double t, tout;
    t = 0.0;
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");
    int Nstep(atoi(argv[1]));
    for (int istep(0); istep < Nstep; ++istep) {
        // output
        out.tabout("# step = ", istep);

        tout = t + dt;

        dvode_(cal_dudt, &NEQ, &u[0], &t, &tout, &itol, &rtol, &atol,
                &itask, &istate, &iopt, &rwork[0], &lrw, &iwork[0], &liw,
                cal_jac, &MF, &rpar[0], &ipar[0]);

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
