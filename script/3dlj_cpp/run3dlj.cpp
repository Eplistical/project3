#include <cmath>
#include <cassert>
#include "lj_particle_3d.hpp"
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
using namespace std;
using ptcl_t = LJ_Particle_3D;

// config
struct Para {
    // basic
    double kT = 1.0;
    double beta = 1.0 / kT;
    double mass = 1.0;
    double mu = -3.2;
    double V = 512.0;
    double L = pow(V, 1.0 / 3.0); 
    double Vc = V;
    double Lc = pow(Vc, 1.0 / 3.0);
    // propagation
    size_t Nstep_eql = 1e5;
    size_t Nstep_run = 7e5;
    // MC related
    double dxmax = 1.0;
    // MD related
    size_t K = 6;
    double dt = 0.005;
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

        vector<double> force(const ptcl_t& p1, const ptcl_t& p2) {
            // return force F12, i.e. force on 1 from 2
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
                    return vector<double>(3, 0.0);
                }
                else {
                    r += dx[i] * dx[i];
                }
            }

            r = sqrt(r);
            if (r > rc) {
                return vector<double>(3, 0.0);
            }
            else {
                double F(LJ_force_raw(r));
                return F / r * dx;
            }
        }

    public:
        double LJ_raw(double r) {
            double tmp(pow(r, -6));
            return 4.0 * tmp * (tmp - 1.0);
        }

        double LJ_force_raw(double r) {
            // return magnitude of force
            double tmp(pow(r, -6));
            return 24.0 / r * tmp * (2.0 * tmp - 1.0);
        }

    public:
        double rc;
        double U_rc;
} LJ;

// arg parser

bool argparse(int argc, char** argv, bool output_flag = true)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("Vc", po::value<double>(&para.Vc), "Vc as a multiple of V")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("K", po::value<double>(&para.K), "MD parameter K")
        ("dt", po::value<double>(&para.dt), "MD time step")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("Nstep_eql", po::value<size_t>(&para.Nstep_eql), "time step for equilibrating the system")
        ("Nstep_run", po::value<size_t>(&para.Nstep_run), "time step for collecting data")
        ;   
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        if (output_flag) {
            ioer::info(desc);
        }
        return false;
    }
    if (vm.count("Vc")) {
        assert(para.Vc <= 1.0 and para.Vc > 0.0);
        para.Vc *= para.V;
        para.Lc = pow(para.Vc, 1.0 / 3.0);
    }
    return true;
}

// helper functions
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


double cal_U(const vector<ptcl_t>& swarm) {
    const size_t N(cal_N(swarm));
    double U(0.0);
    for (size_t i(0); i < N - 1; ++i) {
        for (size_t j(i + 1); j < N; ++j) {
            U += LJ(swarm[i], swarm[j]);
        }
    }
    return U;
}


// main functions
void evolve(vector<ptcl_t>& swarm) {
    // evlove particles with Molecular Dynamics
}


void shuffle(vector<ptcl_t>& swarm) {
    // shuffle particles with Monte Carlo
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
    // create a new particle in Vc
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
    // destcruct an existing particle in Vc
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

void MC_step(size_t istep, vector<ptcl_t>& swarm) {
    // single MC step
    shuffle(swarm);
    if (randomer::rand() < 0.5) {
        create(swarm);
    }
    else {
        destruct(swarm);
    }
}

void MD_step(size_t istep, vector<ptcl_t>& swarm) {
    // a single MD step
    evolve(swarm);
    if (istep % para.K == 0) {
        if (randomer::rand() < 0.5) {
            create(swarm);
        }
        else {
            destruct(swarm);
        }
    }
}


int main(int argc, char** argv) {
    // args
    if (argparse(argc, argv) == false) {
        return 0;
    }
    ioer::info("# V = ", para.V, " Vc = ", para.Vc, " mu = ", para.mu);
    ioer::info("# Nstep_eql = ", para.Nstep_eql, " Nstep_run = ", para.Nstep_run);

    timer::tic();
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
        //ioer::tabout(istep);
        MC_step(istep, swarm);
    }

    // run
    double rho(0.0);
    double U_per_ptcl(0.0);
    for (size_t istep(0); istep < para.Nstep_run; ++istep) {
        rho += cal_N(swarm) / para.V;
        U_per_ptcl += cal_U(swarm) / cal_N(swarm);
        //ioer::tabout(istep, rho / (istep + 1), U_per_ptcl / (istep + 1));
        MC_step(istep, swarm);
    }

    // output
    ioer::tabout(rho / para.Nstep_run, U_per_ptcl / para.Nstep_run);
    ioer::info("# ", timer::toc());

    return 0;
}
