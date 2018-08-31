#ifndef _DVODE_HPP
#define _DVODE_HPP

// c++ interface for DVODE_F90

#include <vector>

extern "C" {
    typedef int int_t;
    typedef void (*cal_dudt_t)(const int_t* /* NEQ */, const double* /* t */, const double* /* u */, double* /* dudt */);
    typedef void (*cal_jac_t)(const int_t* /* NEQ */ , const double* /* t */, const double* /* u */, int_t* /* IA */, int_t* /* JA */, int_t* /* NZ */, double* /* A */);

    extern void dvode_f90_sparse_(  const int_t* /* NEQ */, double* /* u */, double* /* t */, 
                                    double* /* tout */, const int_t* /* itask */, int_t* /* istate */, 
                                    const double* /* rtol */, const double* /* atol */, 
                                    cal_dudt_t /* cal_dudt */, cal_jac_t /* cal_jac */,
                                    const int_t* /* supply_jac */, const int_t* /* reinit */
                                    );
};


inline int_t dvode_sp(std::vector<double>& u, double& t, double& tout, 
        cal_dudt_t cal_dudt, cal_jac_t cal_jac, 
        const double rtol, const double atol, 
        bool reinit = false) 
{
    const int_t N(u.size());
    const int_t reinit_flag(reinit ? 1 : 0);
    const int_t supply_jac_flag(cal_jac ? 1 : 0);
    const int_t itask(1);
    int_t istate(1);

    dvode_f90_sparse_(&N, &u[0], &t, &tout, 
            &itask, &istate, 
            &rtol, &atol, 
            cal_dudt, cal_jac, 
            &supply_jac_flag, &reinit_flag);

    return istate;
}


#endif
