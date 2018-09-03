#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/fmtstring.hpp"
#include "misc/crasher.hpp"
#include "ode.hpp"

using namespace std;

const int Nx(100);
const int Nstep(500);
const double D(1.0);
const double G0(0.0);
const double V0(1.0);
const double V1(-1.0);
const double scan_rate(1.0);
const double dt(abs(V0 - V1) / Nstep / scan_rate);
const double beta(38.92);
const double alpha(0.2);
const double nFA(1.0);
const double kA(0.0);// (0.0);
const double k0(1.0);// (1000.0);

double kf, kb;

const double dx(6.0 * sqrt(2 * abs(V0 - V1) * D / scan_rate) / Nx);
double Dinvdx2(D / dx / dx);
double delta(1.0 / dx);
const double cAinf(1.0);
const double cBinf(0.0);

const int Ntot(Nx * 2);

vector<double> init_u()
{
    vector<double> u(Nx * 2, cBinf);
    for (int i(0); i < Nx; ++i) {
        u[i] = cAinf;
    }
    return u;
}


void cal_dudt(
        const int* /* NEQ */, 
        const double* /* t */, 
        const double* u,
        double* dudt,
        double* /* rpar */,
        int* /* ipar */
        )
{

    // calculate dudt
    dudt[0   ] = Dinvdx2 * (u[1   ] - u[0   ]) - (kf*u[0] - kb*u[Nx]) * delta;// + kA*u[Nx]; 
    dudt[0+Nx] = Dinvdx2 * (u[1+Nx] - u[0+Nx]) + (kf*u[0] - kb*u[Nx]) * delta - kA*u[Nx];

    for (int i(1); i < Nx-1; ++i) {
        dudt[i   ] = Dinvdx2 * (u[i-1   ] - 2*u[i   ] + u[i+1   ]); //+ kA*u[i+Nx];
        dudt[i+Nx] = Dinvdx2 * (u[i-1+Nx] - 2*u[i+Nx] + u[i+1+Nx]) - kA*u[i+Nx];
    }

    dudt[Nx-1   ] = Dinvdx2 * (u[Nx-1   ] - 2*u[Nx-2   ] + u[Nx-3   ]);// + kA*u[Nx-1+Nx];
    dudt[Nx-1+Nx] = Dinvdx2 * (u[Nx-1+Nx] - 2*u[Nx-2+Nx] + u[Nx-3+Nx]) - kA*u[Nx-1+Nx];
}

#if 0
void cal_jac(const int* /* NEQ */, const double* /* t */, const double* u, int* IA, int* JA, int* NZ, double* A)
{
    if (*NZ == 0) {
        // non-zero number
        *NZ = 7 * Nx - 1;
    }
    else {
        // calc IA, JA, A in Fortran's style
        const int ofs(1); // Fortran style
        int count(0);
        IA[0] = ofs;
        for (int j(0); j < Ntot; ++j) {
            if (j == 0) {
                JA[count] = j    + ofs;     A[count] = -Dinvdx2 - kf * delta;   count += 1; // df0/dx0
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                 count += 1; // df1/dx0
                JA[count] = j+Nx + ofs;     A[count] = kf * delta;              count += 1; // df_Nx / dx0
            }
            else if (j < Nx-3) {
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2;              count += 1;
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                 count += 1;
            }
            else if (j == Nx-3) {
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2;              count += 1;
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j+4  + ofs;     A[count] = Dinvdx2;                 count += 1;
            }
            else if (j == Nx-2) {
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2;              count += 1;
                JA[count] = j+1  + ofs;     A[count] = -2*Dinvdx2;              count += 1;
            }
            else if (j == Nx-1) {
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = Dinvdx2;                 count += 1;
            }
            else if (j == Nx) {
                JA[count] = j-Nx + ofs;     A[count] = kb * delta;                      count += 1;
                JA[count] = j    + ofs;     A[count] = -Dinvdx2 - kb * delta - kA;      count += 1;
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                         count += 1;
            }
            else if (j < 2*Nx-3) {
                JA[count] = j-Nx + ofs;     A[count] = kA;                      count += 1;
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2 - kA;         count += 1;
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                 count += 1;
            }
            else if (j == 2*Nx-3) {
                JA[count] = j-Nx + ofs;     A[count] = kA;                      count += 1;
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2 - kA;         count += 1;
                JA[count] = j+1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j+2  + ofs;     A[count] = Dinvdx2;                 count += 1;
            }
            else if (j == 2*Nx-2) {
                JA[count] = j-Nx + ofs;     A[count] = kA;                      count += 1;
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = -2*Dinvdx2 - kA;         count += 1;
                JA[count] = j+1  + ofs;     A[count] = -2*Dinvdx2;              count += 1;
            }
            else if (j == 2*Nx-1) {
                JA[count] = j-Nx + ofs;     A[count] = kA;                      count += 1;
                JA[count] = j-1  + ofs;     A[count] = Dinvdx2;                 count += 1;
                JA[count] = j    + ofs;     A[count] = Dinvdx2 - kA;            count += 1;
            }
            IA[j+1] = count + ofs;
        }
    }
}
#endif

int main(int argc, char** argv) {

    misc::crasher::confirm<>( 2 * D * dt / dx / dx <= 1.0, misc::fmtstring("dt = %f, dx = %f, 2Ddt/dx**2 = %.f > 1.0!", dt, dx, 2*D*dt/dx/dx));

    // output 
    ioer::output_t out("");
    //out.set_precision(10);
    
    // integrator para
    int istate;

    // init
    double V;
    vector<double> u(init_u());
    double t, tout;

    // loop
    V = V0;
    t = 0.0;
    for (int istep(0); istep < Nstep; ++istep) {
        // paras
        tout = t + dt;
        V -= scan_rate * dt;

        kf = k0 * exp(-alpha * beta * (V - G0));
        kb = k0 * exp((1.0 - alpha) * beta * (V - G0));

        // integrate
        /*
        istate = dvode_sp(u, t, tout, cal_dudt, nullptr, 1e-4, 1e-10, true);
        misc::crasher::confirm<>(istate == 2, misc::fmtstring("error: istate = %d", istate));
        */
        dopri5(u, t, tout, cal_dudt, 1e-4, 1e-10);

        // apply boundary condition
        u[Nx-1] = cAinf;
        u[Nx-1+Nx] = cBinf;

        // output
        out.tabout(V, kf*u[0] - kb*u[Nx]);
    }

    return 0;
    /*
    V = V1;
    for (int istep(0); istep < Nstep; ++istep) {
        // paras
        tout = t + dt;
        V += scan_rate * dt;

        kf = k0 * exp(-alpha * beta * (V - G0));
        kb = k0 * exp((1.0 - alpha) * beta * (V - G0));

        // integrate
        istate = dvode_sp(u, t, tout, cal_dudt, nullptr, 1e-4, 1e-10, true);
        misc::crasher::confirm<>(istate == 2, misc::fmtstring("error: istate = %d", istate));

        // apply boundary condition
        u[Nx-1] = cAinf;
        u[Nx-1+Nx] = cBinf;

        // output
        out.tabout(V, kf*u[0] - kb*u[Nx]);
    }
    */

    return 0;
}
