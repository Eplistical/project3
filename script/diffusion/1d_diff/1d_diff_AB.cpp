#include "cvdense.h"
#include "misc/crasher.hpp"
#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include <cmath>
#include <algorithm>

using namespace std;

const int Nx(100);
const size_t Nstep(500);
const double D(1.0);
const double G0(0.0);
const double V0(1.0);
const double V1(-1.0);
const double scan_rate(1.0);
const double dt(abs(V0 - V1) / Nstep / scan_rate);
const double beta(38.92);
const double alpha(0.2);
const double nFA(1.0);
const double kA(0.0);
const double k0(1.0);

double kf(k0 * exp(-alpha * beta * (V0 - G0)));
double kb(k0 * exp((1.0 - alpha) * beta * (V0 - G0)));

const double dx(6.0 * sqrt(2 * abs(V0 - V1) * D / scan_rate) / Nx);

double Dinvdx2(D / dx / dx);
double delta(1.0 / dx);
const double cAinf(1.0);
const double cBinf(0.0);

static void cal_dydt(int /* N */, double /* t */, N_Vector yvec, N_Vector dydtvec, void * /* f_data */) 
{
    const double* y(yvec->data);
    double* dydt(dydtvec->data);

    for (int i(0); i < Nx; ++i) {
        dydt[i] = 0.0; 
        dydt[i+Nx] = 0.0; 
    }

    dydt[0] = Dinvdx2 * (y[1] - y[0]) + kA * y[Nx] - (kf * y[0] - kb * y[Nx]) * delta;
    dydt[Nx] = Dinvdx2 * (y[1+Nx] - y[Nx])- kA * y[Nx] + (kf * y[0] - kb * y[Nx]) * delta;

    for (int i(1); i < Nx-1; ++i) {
        dydt[i] = Dinvdx2 * (y[i-1] - 2 * y[i] + y[i+1]) + kA * y[i+Nx];
        dydt[i+Nx] = Dinvdx2 * (y[i-1+Nx] - 2 * y[i+Nx] + y[i+1+Nx]) - kA * y[i+Nx];
    }

    dydt[Nx-1] = Dinvdx2 * (y[Nx-1] - 2 * y[Nx-2] + y[Nx-3]) + kA * y[Nx-1+Nx];
    dydt[Nx-1+Nx] = Dinvdx2 * (y[Nx-1+Nx] - 2 * y[Nx-2+Nx] + y[Nx-3+Nx]) - kA * y[Nx-1];
}

static void cal_jac(int N, DenseMat jacmat, RhsFn f, void *f_data, double t,
        N_Vector yvec, N_Vector fy, N_Vector ewt, double h, double uround,
        void *jac_data, long int *nfePtr, N_Vector vtemp1,
        N_Vector vtemp2, N_Vector vtemp3) 
{
    const double* y(yvec->data);
    double** jac(jacmat->data); // jac[j][i] = df_i / dx_j, see cvdense.h

    for (int j(0); j < Nx; ++j) {
        for (int i(0); i < Nx; ++i) {
            jac[j][i] = 0.0;
        }
    }

    jac[0][0] = -Dinvdx2 - kf * delta;
    jac[1][0] = Dinvdx2;
    jac[Nx][0] = kA + kb * delta;

    jac[Nx-3][Nx-1] = Dinvdx2;
    jac[Nx-2][Nx-1] = -2 * Dinvdx2;
    jac[Nx-1][Nx-1] = Dinvdx2;
    jac[Nx-1+Nx][Nx-1] = kA;

    jac[Nx][Nx] = -Dinvdx2 - kb * delta - kA;
    jac[Nx+1][Nx] = Dinvdx2;
    jac[0][Nx] = kf * delta;

    jac[Nx-3+Nx][Nx-1+Nx] = Dinvdx2;
    jac[Nx-2+Nx][Nx-1+Nx] = -2 * Dinvdx2;
    jac[Nx-1+Nx][Nx-1+Nx] = Dinvdx2 - kA;

    for (int i(1); i < Nx-1; ++i) {
        jac[i][i] = -2 * Dinvdx2;
        jac[i+1][i] = Dinvdx2;
        jac[i-1][i] = Dinvdx2;
        jac[i+Nx][i] = kA;

        jac[i+Nx][i+Nx] = -2 * Dinvdx2 - kA;
        jac[i+1+Nx][i+Nx] = Dinvdx2;
        jac[i-1+Nx][i+Nx] = Dinvdx2;
    }
}

struct cvode_integrator {
    public:
        int m_NEQ;
        double m_ropt[OPT_SIZE];
        long int m_iopt[OPT_SIZE];
        double m_reltol;
        N_Vector m_abstol;

        double m_t;
        N_Vector m_y;
        void* m_cvode_mem;

    public:
        cvode_integrator(const vector<double>& y, const double T0, const double reltol, const vector<double>& abstol) 
        {
            m_NEQ = y.size();
            m_reltol = reltol;

            m_abstol = N_VNew(m_NEQ, NULL);
            m_y = N_VNew(m_NEQ, NULL);
            m_t = T0;

            copy(abstol.begin(), abstol.begin() + m_NEQ, m_abstol->data);
            copy(y.begin(), y.begin() + m_NEQ, m_y->data);

            m_cvode_mem = CVodeMalloc(m_NEQ, cal_dydt, m_t, m_y, BDF, NEWTON, SV, &m_reltol, m_abstol,
                    NULL, NULL, FALSE, m_iopt, m_ropt, NULL);

            CVDense(m_cvode_mem, cal_jac, NULL);
        }

        cvode_integrator(const vector<double>& y, const double T0, const double reltol, const double abstol) 
        {
            m_NEQ = y.size();
            m_reltol = reltol;

            m_abstol = N_VNew(m_NEQ, NULL);
            m_y = N_VNew(m_NEQ, NULL);
            m_t = T0;

            fill_n(m_abstol->data, m_NEQ, abstol);
            copy(y.begin(), y.begin() + m_NEQ, m_y->data);

            m_cvode_mem = CVodeMalloc(m_NEQ, cal_dydt, m_t, m_y, BDF, NEWTON, SV, &m_reltol, m_abstol,
                    NULL, NULL, FALSE, m_iopt, m_ropt, NULL);

            CVDense(m_cvode_mem, cal_jac, NULL);
        }

        ~cvode_integrator() 
        {
            N_VFree(m_y);
            N_VFree(m_abstol);
            CVodeFree(m_cvode_mem);
        }

    public:
        void integrate(const double T1)
        {
            int flag;

            flag = CVode(m_cvode_mem, T1, m_y, &m_t, NORMAL);
            misc::crasher::confirm<>(flag == SUCCESS, "integrate: flag is not SUCCESS!");
        }

        void reinit(const vector<double>& y, const double T0) 
        {
            m_t = T0;

            copy(y.begin(), y.begin() + m_NEQ, m_y->data);

            CVReInit(m_cvode_mem, cal_dydt, m_t, m_y, BDF, NEWTON, SV, &m_reltol, m_abstol,
                    NULL, NULL, FALSE, m_iopt, m_ropt, NULL);

            CVDense(m_cvode_mem, cal_jac, NULL);
        }

    // observor
    public:
        double get_t() const {
            return m_t;
        }

        vector<double> get_y() const {
            vector<double> y(m_NEQ);
            copy(m_y->data, m_y->data + m_NEQ, y.begin());
            return y;
        }
};

int main(int argc, char** argv) {
    // init
    vector<double> u(Nx * 2, cBinf);
    for (int i(0); i < Nx; ++i) {
        u[i] = cAinf;
    }

    ioer::output_t out("");
    out.set_precision(10);
    cvode_integrator integrator(u, 0.0, 1e-4, 1e-14);

    double V;
    double t;

    // V0 -> V1
    t = 0.0;
    V = V0;
    for (int istep(0); istep < Nstep; ++istep) {
        t += dt;
        V -= scan_rate * dt;

        kf = k0 * exp(-alpha * beta * (V - G0));
        kb = k0 * exp((1.0 - alpha) * beta * (V - G0));

        integrator.integrate(t);

        // apply boundary condition
        u = integrator.get_y();
        u[Nx-1] = cAinf;
        u[Nx-1+Nx] = cBinf;
        integrator.reinit(u, t);
        // output
        out.tabout(V, kf*u[0] - kb*u[Nx]);
    }

    // V1 -> V0
    V = V1;
    for (int istep(0); istep < Nstep; ++istep) {
        t += dt;
        V += scan_rate * dt;

        kf = k0 * exp(-alpha * beta * (V - G0));
        kb = k0 * exp((1.0 - alpha) * beta * (V - G0));

        integrator.integrate(t);

        // apply boundary condition
        u = integrator.get_y();
        u[Nx-1] = cAinf;
        u[Nx-1+Nx] = cBinf;
        integrator.reinit(u, t);
        // output
        out.tabout(V, kf*u[0] - kb*u[Nx]);
    }


    return 0;
}
