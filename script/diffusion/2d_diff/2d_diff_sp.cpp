#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/fmtstring.hpp"
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

const int Nx(100);
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

    using cal_dudt_t = void(*)(int* /* N */, double* /* t */, double* /* u */,
            double* /* dudt */, int* /* ipar */, double* /* rpar */, int* /* ierr */);

    using cal_jac_j_t = void(*)(int* /* N */, double* /* t */, double* /* u */, int* /* j */,
            double* /* jac_j */, int* /* ipar */, double* /* rpar */, int* /* ierr */);

    extern void mebdfso_(
            int* /* N */, double* /* t0 */, double* /* ho */, 
            double* /* u */, double* /* tout */, double* /* tend */,
            int* /* mf */, int* /* idid */, int* /* lout */,
            int* /* lwork */, double* /* work */,
            int* /* liwork */, int* /* iwork */, int* /* maxder */,
            int* /* itol */, double* /* rtol */, double* /* atol */,
            cal_dudt_t /* cal_dudt */, cal_jac_j_t /*cal_jac_j */,
            int* /* nsp */, int* /* nonz */, 
            int* /* ipar */, double* /* rpar */,
            int* /* ierr */
            );
};

void cal_dudt(int* n, double* /* t */, double* u,
        double* dudt, int* /* ipar */, double* /* rpar */, int* /* ierr */)
{
    // calculate dudt
    int idx;

    for (int ix(1); ix < Nx - 1; ++ix) {
        for (int iy(1); iy < Nx - 1; ++iy) {
            idx = ix + Nx * iy;
            dudt[idx] = Dinvdx2 * (u[idx+1] + u[idx-1] + u[idx+Nx] + u[idx-Nx] - 4*u[idx]); 
        }
    }
}

void cal_jac_j(int* /* n */, double* /* t */, double* u, int* j,
        double* jac_j, int* /* ipar */, double* /* rpar */, int* /* ierr */)
{
    // calculate j^th column of jacobian matrix
    int idx(*j);
    int jx, jy;

    jx = idx % Ntot;
    jy = idx / Ntot;

    if (jx != 0 and jx != Nx - 1 and jy != 0 and jy != Nx - 1) {
        jac_j[idx] = -4*Dinvdx2;
        jac_j[idx+1] = Dinvdx2;
        jac_j[idx-1] = Dinvdx2;
        jac_j[idx+Ntot] = Dinvdx2;
        jac_j[idx+Ntot] = Dinvdx2;
    }
}

int main(int argc, char** argv) {
    // output 
    ioer::output_t out("");
    out.set_precision(10);

    // integrator para
    int N(Ntot);

    // tolerance
    int itol(2);
    double rtol(1e-4), atol(1e-14);

    // workspace
    int nonz(Ntot*5), nsp(8*N+2+2*nonz); 
    int liwork(6*N+2*nonz+15), lwork(33*N+2*nonz+nsp+3);
    vector<double> work(lwork);
    vector<int> iwork(liwork);
    iwork[10] = 10000000; // max step

    // parameters pass
    vector<int> ipar(1);
    vector<double> rpar(1);

    // control
    int MF(25); // analytic jacobian
    int maxder(7);
    double ho;
    int idid;
    int lout(0), ierr(0);

    // init
    int Nstep(atoi(argv[1]));
    double t, tout(dt), tend(Nstep * dt);
    t = 0.0;
    vector<double> u = init_u();
    misc::crasher::confirm<>(argc >= 2, "insufficient input para!");

    // loop
    ho = dt;
    idid = 0;
    for (int istep(0); istep < Nstep; ++istep) {
        if (istep == 0) {
            // first call, use idid = 1 to init
            idid = 1;
            mebdfso_(&N, &t, &ho, &u[0], &tout, &tend, &MF, &idid, &lout, 
                    &lwork, &work[0], &liwork, &iwork[0], 
                    &maxder, &itol, &rtol, &atol, 
                    cal_dudt, cal_jac_j, 
                    &nsp, &nonz, &ipar[0], &rpar[0], &ierr);
            idid = 0;
        }

        // output
        out.tabout("# step = ", istep);
        tout = t + dt;

        mebdfso_(&N, &t, &ho, &u[0], &tout, &tend, &MF, &idid, &lout, 
                &lwork, &work[0], &liwork, &iwork[0], 
                &maxder, &itol, &rtol, &atol, 
                cal_dudt, cal_jac_j, 
                &nsp, &nonz, &ipar[0], &rpar[0], &ierr);
        
        misc::crasher::confirm<>(ierr == 0, misc::fmtstring("ierr = %d!!", ierr));

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
