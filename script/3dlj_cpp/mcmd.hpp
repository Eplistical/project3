#ifndef _MCMD_HPP
#define _MCMD_HPP
#include <cmath>
#include <vector>
#include "misc/randomer.hpp"
#include "energy.hpp"

// MC, MD functions

inline bool decide(double x) {
    // generate a prob by e^-x
    if (x > 75.0) 
        return false;
    else if (x < 0.0)
        return true;
    else if (randomer::rand() < exp(-x))
        return true;
    else 
        return false;
}

bool shuffle(std::vector<double>& x, const uint32_t ofs,
        double& U, double& W,
        const double kT, const double dxmax, const double L,
        const double rc, const double Urc)
{
    // shuffle the particle x[ofs:ofs+3] w/ MC algorithm
    uint32_t _3N = x.size();
    assert (ofs % 3 == 0 and ofs < _3N);

    double Uold, Wold, Unew, Wnew;
    double dU, dW;
    std::vector<double> oldx(3);
    one_energy(x, ofs, rc, Urc, L, Uold, Wold);

    for (int k(0); k < 3; ++k) {
        oldx[k] = x[ofs + k];
        x[ofs + k] += randomer::rand(-dxmax, dxmax);
        x[ofs + k] -= L * round(x[ofs + k] / L); // periodic condition
    }
    
    one_energy(x, ofs, rc, Urc, L, Unew, Wnew);

    dU = Unew - Uold;
    dW = Wnew - Wold;

    if (decide(dU / kT)) {
        U += dU;
        W += dW;
        return true;
    }
    else {
        copy(oldx.begin(), oldx.end(), x.begin() + ofs);
        return false;
    }
}

bool create(std::vector<double>& x, const std::vector<double>& newx,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    uint32_t N(_3N / 3);
    double dU, dW;
    double dCB;

    potin(newx, x, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dCB = dU / kT - mu / kT  - log(V / (N + 1));
    if (decide(dCB)) {
        x.resize(_3N + 3);
        copy(newx.begin(), newx.begin() + 3, x.end() - 3);
        U += dU;
        W += dW;
        return true;
    }
    return false;
}

bool destruct(std::vector<double>& x, const uint32_t ofs,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    assert(ofs % 3 == 0 and ofs < _3N);
    uint32_t N(_3N / 3);
    double dU, dW;
    double dDB;

    potout(x, ofs, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dDB = dU / kT + mu / kT - log(N / V);
    if (decide(dDB)) {
        copy(x.end() - 3, x.end(), x.begin() + ofs);
        x.resize(_3N - 3);
        U += dU;
        W += dW;
        return true;
    }
    return false;
}

#endif
