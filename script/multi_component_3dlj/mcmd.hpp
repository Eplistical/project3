#ifndef _MCMD_HPP
#define _MCMD_HPP
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include "misc/randomer.hpp"
#include "energy.hpp"

// MC, MD functions

inline void rand_in_Vc(const uint64_t Ntype, const std::vector<double>& mass, const double kT, const std::vector<double>& Lc, 
        uint64_t& newtype, std::vector<double>& newx, std::vector<double>& newv) 
{
    //newtype = randomer::choice(Ntype);
    newtype = 0;
    //newv = randomer::maxwell_dist(para.mass[newtype], kT);
    newx.resize(3);
    for (int k(0); k < 3; ++k) {
        newx[k] = randomer::rand(-0.5 * Lc[k], 0.5 * Lc[k]);
    }
}

inline bool is_in_Vc(const double* x, const std::vector<double>& Lc) {
    for (int k(0); k < 3; ++k) {
        if (x[k] > 0.5 * Lc[k] or x[k] < -0.5 * Lc[k]) 
            return false;
    }
    return true;
}

inline std::vector<uint64_t> get_Vc_idx(const std::vector<double>& x, const std::vector<double>& Lc) 
{
    // get indices for particles that are in the control colume
    const uint64_t Ntot(x.size() / 3);

    std::vector<uint64_t> Vc_idx;
    Vc_idx.reserve(Ntot);

    for (uint64_t i(0); i < Ntot; ++i) {
        if (is_in_Vc(&x[3 * i], Lc)) {
            Vc_idx.push_back(i);
        }
    }
    return Vc_idx;
}

inline uint64_t get_Nc(const uint64_t thetype, const std::vector<uint64_t>& type, 
        const std::vector<double>& x, const std::vector<double>& Lc)
{
    // get ptcl number for thetype molecule in the control volume
    uint64_t Nc(0);
    const uint64_t Ntot(type.size());
    for (uint64_t i(0); i < Ntot; ++i) {
        if ((type[i] == thetype) and (is_in_Vc(&x[3 * i], Lc))) {
            Nc += 1;
        }
    }
    return Nc;
}

inline std::vector<uint64_t> get_N(const uint64_t Ntype, const std::vector<uint64_t>& type)
{
    // get ptcl number for all types of molecules
    std::vector<uint64_t> N(Ntype);
    for (const uint64_t& i : type) {
        N[i] += 1;
    }
    return N;
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

bool shuffle(std::vector<double>& x, const uint64_t idx,
        const uint64_t Ntype, const std::vector<uint64_t> type, 
        const double kT, const double dxmax,
        const std::vector<double>& sigma, const std::vector<double>& epsilon, 
        const std::vector<double>& rc, const std::vector<double>& Urc, 
        const std::vector<double>& L, 
        const std::vector<double>& ULRC, const std::vector<double>& WLRC,
        double& U, double& W)
{
    // shuffle the particle idx w/ MC algorithm
    double U0, W0, U1, W1;
    double dU, dW;

    std::vector<double> oldx(3);
    one_energy(x, idx, Ntype, type, sigma, epsilon, rc, Urc, L, U0, W0, nullptr);

    for (int k(0); k < 3; ++k) {
        oldx[k] = x[k + 3 * idx];
        x[k + 3 * idx] += randomer::rand(-dxmax, dxmax);
        x[k + 3 * idx] -= L[k] * round(x[k + 3 * idx] / L[k]); // periodic condition
    }
    
    one_energy(x, idx, Ntype, type, sigma, epsilon, rc, Urc, L, U1, W1, nullptr);

    dU = U1 - U0;
    dW = W1 - W0;

    if (decide(dU / kT)) {
        U += dU;
        W += dW;
        return true;
    }
    else {
        copy(oldx.begin(), oldx.end(), x.begin() + 3 * idx);
        return false;
    }
}

bool create(std::vector<double>& x, std::vector<double>& v, 
        std::vector<uint64_t>& type,
        const std::vector<double>& newx, const std::vector<double>& newv, 
        const uint64_t newtype,
        const uint64_t Ntype, const vector<uint64_t>& N,
        const std::vector<double>& sigma, const std::vector<double>& epsilon, 
        const std::vector<double>& rc, const std::vector<double>& Urc,
        const std::vector<double>& L, const double kT, const std::vector<double>& mu, 
        const double Vc, const uint64_t Nc,
        const std::vector<double>& ULRC, const std::vector<double>& WLRC,
        double& U, double& W)
{
    // attempt to create a new particle at newx w/ velocity newv
    const uint64_t Ntot(x.size() / 3);
    double dU, dW;
    double dCB;

    potin(x, newx, type, newtype, Ntype, N,
            sigma, epsilon, rc, Urc, L, ULRC, WLRC, dU, dW);

    dCB = dU / kT - mu[newtype] / kT  - log(Vc / (Nc + 1));
    if (decide(dCB)) {
        type.push_back(newtype);
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
        std::vector<uint64_t>& type,
        const uint64_t idx,
        const uint64_t Ntype, const vector<uint64_t>& N,
        const std::vector<double>& sigma, const std::vector<double>& epsilon, 
        const std::vector<double>& rc, const std::vector<double>& Urc,
        const std::vector<double>& L, const double kT, const std::vector<double>& mu, 
        const double Vc, const uint64_t Nc,
        const std::vector<double>& ULRC, const std::vector<double>& WLRC,
        double& U, double& W)
{
    // attempt to destruct an existing particle idx
    const uint64_t idxtype(type[idx]);
    double dU, dW;
    double dDB;

    potout(x, idx, type, Ntype, N,
            sigma, epsilon, rc, Urc, L, ULRC, WLRC, dU, dW);

    dDB = dU / kT + mu[idxtype] / kT - log(Nc / Vc);
    if (decide(dDB)) {
        type.erase(type.begin() + idx);
        x.erase(x.begin() + idx * 3, x.begin() + idx * 3 + 3);
        U += dU;
        W += dW;

        if (not v.empty()) {
            v.erase(v.begin() + idx * 3, v.begin() + idx * 3 + 3);
        }
        return true;
    }
    return false;
}


/*
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

    // velocity verlet step 1
    F.assign(_3N, 0.0);
    all_energy(x, rc, Urc, L, ULRC0, WLRC0, U, W, &F[0], true);
    x = x + v * dt + 0.5 / mass * dt * dt * F;
    v = v + 0.5 / mass * dt * F;

    // periodic condition
    for (auto& xi : x) {
        xi -= L * round(xi / L);
    }

    // velocity verlet step 2
    F.assign(_3N, 0.0);
    all_energy(x, rc, Urc, L, ULRC0, WLRC0, U, W, &F[0], true);
    v = v + 0.5 / mass * dt * F;

    for (auto& xi : x) {
        if (xi > L / 2 or xi < -L / 2) {
            ioer::info("out of range!!!!");
        }
    }
}
*/

#endif
