#include <cmath>
#include <cassert>
#include "lj_particle_3d.hpp"
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
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


size_t cal_N(const vector<ptcl_t>& swarm) {
    return swarm.size();
}

bool is_in_Vc(const ptcl_t& ptcl) {
    for (double xi : ptcl.x) {
        if (xi > para.Lc or xi < 0.0) {
            return false;
        }
    }
    return true;
}

vector<size_t> get_idx_in_Vc(const vector<ptcl_t>& swarm) {
    const size_t N(cal_N(swarm));
    vector<size_t> idx_in_Vc;
    idx_in_Vc.reserve(N);
    for (size_t i(0); i < N; ++i) {
        if (is_in_Vc(swarm[i])) {
            idx_in_Vc.push_back(i);
        }
    }
    return idx_in_Vc;
}

size_t cal_Nc(const vector<ptcl_t>& swarm) {
    const size_t N(cal_N(swarm));
    size_t Nc(0);
    for (size_t i(0); i < N; ++i) {
        if (is_in_Vc(swarm[i])) {
            ++Nc;
        }
    }
    return Nc;
}

ptcl_t rand_in_Vc() {
    return ptcl_t(
            randomer::vrand(3, 0.0, para.Lc),
            randomer::maxwell_dist<>(para.mass, para.kT),
            para.mass,
            para.kT
            );
}

void shuffle(vector<ptcl_t>& swarm) {
    const size_t N(cal_N(swarm));
    for (size_t i(0); i < N; ++i) {
        vector<double> dx = randomer::vrand(3, -para.dxmax, para.dxmax);
        ptcl_t last_ptcl(swarm[i]);

        // move w/ periodic condition
        swarm[i].x = swarm[i].x + dx;
        for (int k(0); k < 3; ++k) {
            if (swarm[i].x[k] < 0.0) {
                swarm[i].x[k] += para.L;
            }
            else if (swarm[i].x[k] > para.L) {
                swarm[i].x[k] -= para.L;
            }
            assert(swarm[i].x[k] >= 0.0 and swarm[i].x[k] <= para.L);
        }
        // calc U
        double dU(0.0);
        for (size_t j(0); j < N; ++j) {
            if (i != j) {
                dU += LJ(swarm[j], swarm[i]) - LJ(swarm[j], last_ptcl);
            }
        }
        // make decision
        const double prob(exp(-para.beta * dU));
        if (randomer::rand() >= prob) {
            swarm[i].x = last_ptcl.x;
        }
    }
}

void create(vector<ptcl_t>& swarm) {
    ptcl_t new_ptcl(rand_in_Vc());
    // calc dU
    const size_t N(cal_N(swarm));
    double dU(0.0);
    for (int i(0); i < N; ++i) {
        dU += LJ(swarm[i], new_ptcl);
    }
    // make decision
    const size_t Nc(cal_Nc(swarm));
    const double prob(exp(-para.beta * (dU - para.mu)) * para.Vc / (Nc + 1));
    if (randomer::rand() < prob) {
        swarm.push_back(new_ptcl);
    }
}

void destruct(vector<ptcl_t>& swarm) {
    vector<size_t> idx_in_Vc(get_idx_in_Vc(swarm));
    if (not idx_in_Vc.empty()) {
        const size_t Nc(idx_in_Vc.size());
        const size_t idx(randomer::choice(idx_in_Vc));
        // calc dU
        const size_t N(cal_N(swarm));
        double dU(0.0);
        for (int i(0); i < N; ++i) {
            if (idx != i) {
                dU -= LJ(swarm[i], swarm[idx]);
            }
        }
        // make decision
        const double prob(exp(-para.beta * (dU + para.mu)) * Nc / para.Vc);
        if (randomer::rand() < prob) {
            swarm.erase(swarm.begin() + idx);
        }
    }
}

void MC_step(vector<ptcl_t>& swarm) {
    shuffle(swarm);
    if (randomer::rand() < 0.5) {
        create(swarm);
    }
    else {
        destruct(swarm);
    }
}



int main(int argc, char** argv) {
    // init
    vector<ptcl_t> swarm {
        ptcl_t(
                vector<double>(3, 0.5 * para.L),
                randomer::maxwell_dist(para.mass, para.kT),
                para.mass,
                para.kT
              )
    };

    // equilibrate
    for (size_t istep(0); istep < para.Nstep_eql; ++istep) {
        ioer::tabout(istep);
        MC_step(swarm);
    }

    // run
    double rho(0.0);
    for (size_t istep(0); istep < para.Nstep_run; ++istep) {
        rho += cal_N(swarm) / para.V;
        ioer::tabout(istep, rho / (istep + 1));
        MC_step(swarm);
    }

    return 0;
}
