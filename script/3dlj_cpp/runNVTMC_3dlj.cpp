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
    const size_t N = 108;
    const double V = N / 0.9;
    const double L = pow(V, 1.0 / 3.0); 
    const double rc = 2.5;

    double kT = 0.9;
    double beta = 1.0 / kT;
    double mass = 1.0;
    // propagation
    size_t Nstep_eql = 1e5;
    size_t Nstep_run = 7e5;
    // MC related
    double dxmax = 1.0;
} para;

// potential
struct LJ_shifted {
    public:
        LJ_shifted(double RC, double LBOX) noexcept 
            : rc(RC), Lbox(LBOX)
        {
            // -- basic parameters -- //
            rc2 = rc * rc;
            Vbox = Lbox * Lbox * Lbox;

            // -- shifted energy -- //
            U_rc = 4.0 * (pow(rc, -12) - pow(rc, -6));

            
            // -- long range correction -- //
            
            // U_LRC0 = U_LRC / N**2
            U_LRC0 = 8.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 3.0 * pow(rc, -3)); 
            // P_LRC0 = P_LRC / N**2
            P_LRC0 = 32.0 / 9.0 * M_PI / Vbox / Vbox * (pow(rc, -9) - 1.5 * pow(rc, -3));
            // mu_LRC0 = mu_LRC / N
            mu_LRC0 = 16.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 3.0 * pow(rc, -3));
        }

    public:
        double operator()(const vector<double>& x1, const vector<double>& x2) {
            /**
             * calculate shifted potential, virial, and force between coordinates 1 & 2
             * return shifted potential
             * store virial and force in public variables
             */
            vector<double> dx(x1 - x2);
            double r2(0.0);
            for (auto& dxi : dx) {
                // min imag
                dxi -= Lbox * round(dxi / Lbox);
                assert(dxi <= 0.5 * Lbox and dxi >= -0.5 * Lbox);

                r2 += dxi * dxi;
                if (r2 > rc2) {
                    U_shifted = 0.0;
                    virial = 0.0;
                    force.assign(3, 0.0);
                    return U_shifted;
                }
            }

            // calc U_shifted, virial & force
            const double r6(r2 * r2 * r2);
            const double r6_inv(1.0 / r6);
            U_shifted = 4.0 * r6_inv * (r6_inv - 1.0) - U_rc;
            virial = 24.0 * r6_inv * (2.0 * r6_inv - 1.0);
            force = virial / r2 * dx;

            return U_shifted;
        }

    public:
        double U_shifted;
        double virial;
        vector<double> force;

    public:
        double rc, rc2, U_rc;
        double Lbox, Vbox;
        double U_LRC0, P_LRC0, mu_LRC0;
} LJ(para.rc, para.L);

// arg parser

bool argparse(int argc, char** argv, bool output_flag = true)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
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
    return true;
}

// helper functions
void boundary_condition(vector<double>& x) {
    for (auto& xi : x) {
        if (xi < 0.0) {
            xi += para.L;
        }
        else if (xi > para.L) {
            xi -= para.L;
        }
        assert(xi >= 0.0 and xi <= para.L);
    }
}

double cal_U(const vector<ptcl_t>& swarm) {
    const size_t N(swarm.size());
    double U(0.0);
    for (size_t i(0); i < N - 1; ++i) {
        for (size_t j(i + 1); j < N; ++j) {
            U += LJ(swarm[i].x, swarm[j].x);
        }
    }
    return U;
}

double cal_W(const vector<ptcl_t>& swarm) {
    const size_t N(swarm.size());
    double W(0.0);
    for (size_t i(0); i < N - 1; ++i) {
        for (size_t j(i + 1); j < N; ++j) {
            LJ(swarm[i].x, swarm[j].x);
            W += LJ.virial / 3.0;
        }
    }
    return W;
}

// main functions

size_t shuffle(vector<ptcl_t>& swarm) {
    // shuffle particles with Monte Carlo
    size_t Naccept(0);
    const size_t N(swarm.size());
    for (size_t i(0); i < N; ++i) {
        vector<double> dx = randomer::vrand(3, -para.dxmax, para.dxmax);
        ptcl_t last_ptcl(swarm[i]);
        // move w/ boundary condition
        swarm[i].x = swarm[i].x + dx;
        boundary_condition(swarm[i].x);
        // calc U
        double dU(0.0);
        for (size_t j(0); j < N; ++j) {
            if (i != j) {
                dU += LJ(swarm[j].x, swarm[i].x) - LJ(swarm[j].x, last_ptcl.x);
            }
        }
        // make decision
        const double prob(exp(-para.beta * dU));
        if (randomer::rand() >= prob) {
            swarm[i].x = last_ptcl.x;
        }
        else {
            Naccept += 1;
        }
    }
    return Naccept;
}

vector<ptcl_t> init_swarm(int N) {
    vector<ptcl_t> swarm;
    swarm.reserve(N);
    for (int i(0); i < N; ++i) {
        swarm.push_back(
                ptcl_t(
                    randomer::vrand(3, 0.0, para.L),
                    randomer::maxwell_dist(para.mass, para.kT),
                    para.mass,
                    para.kT
                    )
                );
    }
    return swarm;
}

int main(int argc, char** argv) {
    // args
    if (argparse(argc, argv) == false) {
        return 0;
    }
    timer::tic();
    // setup
    vector<ptcl_t> swarm(init_swarm(para.N));

    // equilibrate
    for (size_t istep(1); istep <= para.Nstep_eql; ++istep) {
        ioer::info(istep);
        shuffle(swarm);
    }

    // sampling
    double Utot(0.0);
    double P(0.0);
    size_t Naccept(0);
    for (size_t istep(1); istep <= para.Nstep_run; ++istep) {
        size_t N(swarm.size());
        Utot += cal_U(swarm) + LJ.U_LRC0 * N * N;
        P += N / para.V * para.kT + cal_W(swarm) / para.V + LJ.P_LRC0 * N * N;
        ioer::tabout(istep, Utot / istep / N, P / istep);
        Naccept += shuffle(swarm);
    }

    // output
    ioer::info("# rc = ", para.rc, " kT = ", para.kT);
    ioer::info("# V = ", para.V, " rho = ", para.N / para.V, " N = ", para.N);
    ioer::info("# Nstep_eql = ", para.Nstep_eql, " Nstep_run = ", para.Nstep_run);

    ioer::info("# acc ratio = ", static_cast<double>(Naccept) / para.Nstep_run / swarm.size());
    ioer::info("# U = ", Utot / para.Nstep_run / para.N);
    ioer::info("# P = ", P / para.Nstep_run);
    ioer::info("# acc ratio = ", static_cast<double>(Naccept) / swarm.size() / para.Nstep_run);
    ioer::info("# ", timer::toc());

    return 0;
}
