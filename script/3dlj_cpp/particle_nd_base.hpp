#ifndef _PARTICLE_ND_BASE_HPP
#define _PARTICLE_ND_BASE_HPP
// module for N-dimensional particle base class (abstract class)

#include <string>
#include <vector>
#include <numeric>

namespace {
    struct Particle_ND
    {
        public:
            Particle_ND(   const std::vector<double>& X, const std::vector<double>& V, double MASS, double KT) noexcept :
                x(X), v(V), mass(MASS), kT(KT), mass_inv(1.0 / MASS), kT_inv(1.0 / KT), dim(X.size())
                {
                    v.resize(dim);
                }

            virtual ~Particle_ND() noexcept = default;

        public:
            double get_Ek() const noexcept { 
                return 0.5 * mass * std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
            }

        public:
            void evolve(double dt, const std::string& alg = "verlet") {
                do_evolve(dt, alg); 
            }

        public:
            size_t dim;
            double mass, mass_inv;
            double kT, kT_inv;
            std::vector<double> x, v;

        private:
            virtual void do_evolve(double dt, const std::string& alg) = 0;
    };
};

#endif // _PARTICLE_ND_BASE_HPP
