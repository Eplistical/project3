#include <cmath>
#include <cassert>
#include "lj_particle_3d.hpp"
#include "misc/vector.hpp"
#include "misc/ioer.hpp"

using namespace std;
using ptcl_t = LJ_Particle_3D;

// config
struct Para {
    double kT = 1.0;
    double beta = 1.0 / kT;
    double mass = 1.0;
    double mu = -3.2;
    double V = 512.0;
    double L = pow(V, 1.0 / 3.0); 
    double Vc = V;
    double Lc = pow(Vc, 1.0 / 3.0);

    double dxmax = 1.0;
    int Nstep_eql = 1e5;
    int Nstep_run = 3e5;
} para;

// potential
struct LJ_with_cutoff {
    public:
        LJ_with_cutoff() noexcept {
            rc = 2.5;
            U_rc = LJ_raw(rc);
        }

    public:
        double operator()(const ptcl_t& p1, const ptcl_t& p2) {
            vector<double> dx(p1.x - p2.x);
            double r(0.0);
            for (int i(0); i < 3; ++i) {
                // minimum image
                if (dx[i] > 0.5 * para.L) {
                    dx[i] = para.L - dx[i];
                } 
                else if (dx[i] < -0.5 * para.L) {
                    dx[i] = para.L + dx[i];
                }
                assert(dx[i] <= 0.5 * para.L and dx[i] >= -0.5 * para.L);

                if (abs(dx[i]) > rc) {
                    return 0.0;
                }
                else {
                    r += dx[i] * dx[i];
                }
            }

            r = sqrt(r);
            if (r > rc) {
                return 0.0;
            }
            else {
                return LJ_raw(r) - U_rc;
            }
        }

    public:
        double LJ_raw(double r) {
            double tmp(pow(r, -6));
            return 4.0 * tmp * (tmp - 1.0);
        }

    public:
        double rc;
        double U_rc;
} LJ;



int main(int argc, char** argv) {
    vector<double> zzz{0,0,0};
    ptcl_t p1(zzz, zzz, 1.0, 1.0);
    ptcl_t p2(zzz, zzz, 1.0, 1.0);

    for (double x(0.8); x < 5.0; x += 0.02)
    {
        p2.x[1] = x;
        ioer::tabout(x, LJ(p1, p2), LJ.LJ_raw(x), LJ.U_rc);
    }

    return 0;
}
