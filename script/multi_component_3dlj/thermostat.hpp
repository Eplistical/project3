#ifndef _THERMOSTAT_HPP
#define _THERMOSTAT_HPP
#include <cmath>
#include <algorithm>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"


void berendsen_thermostat(std::vector<double>& v,
        const double mass, const double kT,
        const double tau, const double dt) 
{
    const double kTnow(mass * mean(v * v));
    const double lambda(sqrt(1 + dt / tau * (kT / kTnow - 1)));
    v = v * lambda;
}


void andersen_thermostat(std::vector<double>& v,
        const uint64_t Ntype, const std::vector<uint64_t>& type,
        const std::vector<double>& mass, const double kT,
        const double nu, const double dt) 
{
    const double nudt(nu * dt);
    const uint64_t Ntot(type.size());
    std::vector<double> vnew;

    for (uint64_t i(0); i < Ntot; ++i) {
        if (randomer::rand() < nudt) {
            vnew = randomer::maxwell_dist(mass[type[i]], kT);
            std::copy(vnew.begin(), vnew.begin() + 3, v.begin() + i * 3);
        }
    }
}


#endif
