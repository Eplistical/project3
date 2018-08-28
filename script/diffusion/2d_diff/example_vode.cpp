#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "misc/crasher.hpp"
#include "misc/fmtstring.hpp"
using namespace std;

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
    dudt[0] = -0.04*u[0] + 1e4*u[1]*u[2];
    dudt[1] = 0.04*u[0] - 1e4*u[1]*u[2] - 3e7*u[1]*u[1];
    dudt[2] = 3e7*u[1]*u[1];
}

void cal_jac(int* NEQ, double* /* t */, 
        double* u, 
        int* /* ML */, int* /* MU */, 
        double* jac,  int* /* NRPD */,
        double* /* rpar */, int* /* ipar */) 
{
    int N(*NEQ);
    jac[0 + 0*N] = -0.04;
    jac[0 + 1*N] = 1e4*u[2];
    jac[0 + 2*N] = 1e4*u[1];
    jac[1 + 0*N] = 0.04;
    jac[1 + 1*N] = -1e4*u[2] - 6e7*u[1];
    jac[1 + 2*N] = -1e4*u[1];
    jac[2 + 0*N] = 0.0;
    jac[2 + 1*N] = 6e7*u[1];
    jac[2 + 2*N] = 0.0;
}

int main(int argc, char** argv) 
{
    int NEQ(3);
    int MF(21); // BDF
    double rtol(1e-4);
    double atol(1e-14);
    int itol(1), itask(1), iopt(0);
    int lrw(22 + 9*NEQ + 2*NEQ*NEQ);
    int liw(30 + NEQ);
    vector<double> rwork(lrw, 0.0);
    vector<int> iwork(liw, 0.0);
    vector<double> rpar(10000, 0.0);
    vector<int> ipar(10000, 0);
    int istate(1);

    vector<double> u { 1.0, 0.0, 0.0 };
    double t(0.0), tout(0.4);

    for (int it(0); it < 11; ++it) {
        dvode_(cal_dudt, &NEQ, &u[0], &t, &tout, &itol, &rtol, &atol,
                &itask, &istate, &iopt, &rwork[0], &lrw, &iwork[0], &liw,
                cal_jac, &MF, &rpar[0], &ipar[0]);
        misc::crasher::confirm<>(istate == 2, misc::fmtstring("error: istate = %d", istate));
        ioer::tabout(misc::fmtstring("t = %12.4e  u = %15.6e%15.6e%15.6e", tout, u[0], u[1], u[2]));
        tout *= 10.0;
    }
    return 0;
}
