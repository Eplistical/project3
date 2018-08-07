#ifndef _MCMD_HPP
#define _MCMD_HPP
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include "misc/randomer.hpp"
#include "energy.hpp"

// MC, MD functions

inline bool is_in_Vc(const double* x, const double Lc) {
    // check if x[0:3] is in the control volume [-0.5*Lc, 0.5*Lc]
    for (int k(0); k < 3; ++k) {
        if (x[k] > 0.5 * Lc or x[k] < -0.5 * Lc) {
            return false;
        }
    }
    return true;
}

inline std::vector<uint64_t> get_Vc_idx(const std::vector<double>& x, const double Lc) {
    // get indices for particles that are in the control colume
    const uint64_t _3N(x.size());
    const uint64_t N(_3N / 3);

    std::vector<uint64_t> Vc_idx;
    Vc_idx.reserve(N);
    for (uint64_t i(0); i < N; ++i) {
        if (is_in_Vc(&x[3 * i], Lc)) {
            Vc_idx.push_back(i);
        }
    }
    return Vc_idx;
}

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

bool shuffle(std::vector<double>& x, const uint64_t ofs,
        double& U, double& W,
        const double kT, const double dxmax, const double L,
        const double rc, const double Urc,
        const double ULRC0, const double WLRC0)
{
    // shuffle the particle x[ofs:ofs+3] w/ MC algorithm
    const uint64_t _3N = x.size();
    assert (ofs % 3 == 0 and ofs < _3N);

    double U0, W0, U1, W1;
    double dU, dW;
    std::vector<double> F;

    std::vector<double> oldx(3);
    one_energy(x, ofs, rc, Urc, L, U0, W0, &F[0]);

    for (int k(0); k < 3; ++k) {
        oldx[k] = x[ofs + k];
        x[ofs + k] += randomer::rand(-dxmax, dxmax);
        x[ofs + k] -= L * round(x[ofs + k] / L); // periodic condition
    }
    
    one_energy(x, ofs, rc, Urc, L, U1, W1, &F[0]);

    dU = U1 - U0;
    dW = W1 - W0;

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

bool create(std::vector<double>& x, std::vector<double>& v, 
        const std::vector<double>& newx, const std::vector<double>& newv,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    // attempt to create a new particle at newx w/ velocity newv
    const uint64_t _3N(x.size());
    const uint64_t N(_3N / 3);
    double dU, dW;
    double dCB;

    potin(x, newx, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dCB = dU / kT - mu / kT  - log(V / (N + 1));
    if (decide(dCB)) {
        x.insert(x.end(), newx.begin(), newx.begin() + 3);
        U += dU;
        W += dW;
        if (not v.empty()) {
            v.insert(v.end(), newv.begin(), newv.begin() + 3);
        }
        return true;
    }
    return false;
}

bool destruct(std::vector<double>& x, std::vector<double>& v,
        const uint64_t ofs,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    // attempt to destruct an existing particle x[ofs]
    // if success, rescale rest v to maintain total momentum
    const uint64_t _3N(x.size());
    const uint64_t N(_3N / 3);
    assert(ofs % 3 == 0 and ofs < _3N);
    double dU, dW;
    double dDB;

    potout(x, ofs, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dDB = dU / kT + mu / kT - log(N / V);
    if (decide(dDB)) {
        x.erase(x.begin() + ofs, x.begin() + ofs + 3);
        U += dU;
        W += dW;

        if (not v.empty()) {
            /*
            // calc coefficient
            std::vector<double> coef(3, 0.0);
            for (uint64_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
                if (ofs != ofs_i) {
                    for (int k(0); k < 3; ++k) {
                        coef[k] += std::abs(v[ofs_i + k]);
                    }
                }
            }
            for (int k(0); k < 3; ++k) {
                coef[k] = v[ofs + k] / coef[k];
            }
            // rescale rest v
            for (uint64_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
                if (ofs != ofs_i) {
                    for (int k(0); k < 3; ++k) {
                        v[ofs_i + k] += std::abs(v[ofs_i + k]) * coef[k];
                    }
                }
            }
            // erase v[ofs]
            v.erase(v.begin() + ofs, v.begin() + ofs + 3);
            */
            v.erase(v.begin() + ofs, v.begin() + ofs + 3);
            //v = v * sqrt(mean(v * v) / kT);
        }
        return true;
    }
    return false;
}


void evolve(std::vector<double>& x, std::vector<double>& v,
        double& U, double& W, const double L,
        const double dt, const double mass, 
        const double kT, const double nu,
        const double rc, const double Urc,
        const double ULRC0, const double WLRC0)
{
    // evolve the system w/ MD algorithm
    static std::vector<double> F;
    const uint64_t _3N(x.size());
    const uint64_t _3N_3(_3N - 3);
    const uint64_t N(_3N / 3);
    double Uij, Wij;
    F.resize(_3N);

    // velocity verlet
    F.assign(_3N, 0.0);
    all_energy(x, rc, Urc, L, ULRC0, WLRC0, U, W, &F[0], true);
    x = x + v * dt + 0.5 / mass * dt * dt * F;
    v = v + 0.5 / mass * dt * F;

    F.assign(_3N, 0.0);
    all_energy(x, rc, Urc, L, ULRC0, WLRC0, U, W, &F[0], true);
    v = v + 0.5 / mass * dt * F;

    // Andersen Thermostat
    const double nudt(nu * dt);
    std::vector<double> vnew;
    for (uint64_t ofs(0); ofs < _3N; ofs += 3) {
        if (randomer::rand() < nudt) {
            vnew = randomer::maxwell_dist(mass, kT);
            std::copy(vnew.begin(), vnew.begin() + 3, v.begin() + ofs);
        }
    }
}

#endif
