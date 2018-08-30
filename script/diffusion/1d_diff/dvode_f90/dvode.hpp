#ifndef _DVODE_HPP
#define _DVODE_HPP

// c++ interface for DVODE_F90

#include <vector>

extern "C" {
    typedef void (*cal_dudt_t)(const int* /* NEQ */, const double* /* t */, const double* /* u */, double* /* dudt */);
    typedef void (*cal_jac_t)(const int* /* NEQ */ , const double* /* t */, const double* /* u */, int* /* IA */, int* /* JA */, int* /* NZ */, double* /* A */);

    extern void dvode_f90_sparse_(  const int* /* NEQ */, double* /* u */, double* /* t */, 
                                    double* /* tout */, const int* /* itask */, int* /* istate */, 
                                    const double* /* rtol */, const double* /* atol */, 
                                    cal_dudt_t /* cal_dudt */, cal_jac_t /* cal_jac */,
                                    const int* /* supply_jac */, const int* /* reinit */
                                    );
};


inline int dvode_sp(std::vector<double>& u, double& t, double& tout, 
        cal_dudt_t cal_dudt, cal_jac_t cal_jac, 
        const double rtol, const double atol, 
        bool reinit = false) 
{
    const int N(u.size());
    const int reinit_flag(reinit ? 1 : 0);
    const int supply_jac_flag(cal_jac ? 1 : 0);
    const int itask(1);
    int istate(1);

    dvode_f90_sparse_(&N, &u[0], &t, &tout, 
            &itask, &istate, 
            &rtol, &atol, 
            cal_dudt, cal_jac, 
            &supply_jac_flag, &reinit_flag);

    return istate;
}


#endif
