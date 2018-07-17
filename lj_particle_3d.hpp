#ifndef _LJ_PARTICLE_3D_HPP
#define _LJ_PARTICLE_3D_HPP

#include <string>
#include <cassert>
#include "particle_nd_base.hpp"

namespace {
    struct LJ_Particle_3D final : public Particle_ND {
        public:
            LJ_Particle_3D(const std::vector<double>& X, const std::vector<double>& V, double MASS, double KT) 
                : Particle_ND(X, V, MASS, KT)
            {
                assert(this->dim == 3);
            }

            ~LJ_Particle_3D() noexcept = default;

        private:
            void do_evolve(double dt, const std::string& alg) override
            {
            }
    };
};

#endif // _LJ_PARTICLE_3D_HPP
