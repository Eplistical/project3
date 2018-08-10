#ifndef _THERMOSTAT_HPP
#define _THERMOSTAT_HPP
#include <cmath>
#include <cassert>
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
        const double mass, const double kT,
        const double nu, const double dt) 
{
    const double nudt(nu * dt);
    const uint64_t _3N(v.size());
    std::vector<double> vnew;

    for (uint64_t ofs(0); ofs < _3N; ofs += 3) {
        if (randomer::rand() < nudt) {
            vnew = randomer::maxwell_dist(mass, kT);
            std::copy(vnew.begin(), vnew.begin() + 3, v.begin() + ofs);
        }
    }
}


#endif
