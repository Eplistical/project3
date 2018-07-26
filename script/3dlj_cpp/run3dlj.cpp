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
    const double V = 120.0;
    const double L = pow(V, 1.0 / 3.0); 
    const double rc = 3.0;

    double kT = 0.9;
    double beta = 1.0 / kT;
    double mass = 1.0;
    double mu = -2.5;

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
        double operator()(const ptcl_t& p1, const ptcl_t& p2) {
            /**
             * calculate shifted potential, virial, and force between particles 1 & 2
             * return shifted potential
             * store virial and force in public variables
             */
            vector<double> dx(p1.x - p2.x);
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
        ("Vc", po::value<double>(&para.Vc), "Vc as a multiple of V")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("K", po::value<size_t>(&para.K), "MD parameter K")
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

void boundary_condition(vector<double>& x) {
    for (int k(0); k < 3; ++k) {
        if (x[k] < 0.0) {
            x[k] += para.L;
        }
        else if (x[k] > para.L) {
            x[k] -= para.L;
        }
        assert(x[k] >= 0.0 and x[k] <= para.L);
    }
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

vector< vector<double> > cal_force(const vector<ptcl_t>& swarm) {
    const size_t N(cal_N(swarm));
    vector< vector<double> > force(N, vector<double>(3, 0.0));
    for (size_t i(0); i < N - 1; ++i) {
        for (size_t j(i + 1); j < N; ++j) {
            LJ(swarm[i], swarm[j]);
            force[i] = force[i] + LJ.force;
            force[j] = force[j] - LJ.force;
        }
    }
    return force;
}

// main functions
void evolve(vector<ptcl_t>& swarm) {
    // evlove particles with Molecular Dynamics
    const size_t N(cal_N(swarm));
    vector< vector<double> > force;

    // 
    force = cal_force(swarm);
    for (size_t i(0); i < N; ++i) {
        // adjust v
        swarm[i].v = swarm[i].v + para.dt * swarm[i].mass_inv * force[i];
        // move w/ boundary condition
        swarm[i].x = swarm[i].x + para.dt * swarm[i].v;
        boundary_condition(swarm[i].x);
    }
    // velocity verlet
    force = cal_force(swarm);
    for (size_t i(0); i < N; ++i) {
        swarm[i].v = swarm[i].v + 0.5 * para.dt * swarm[i].mass_inv * force[i];
        // move w/ boundary condition
        swarm[i].x = swarm[i].x + para.dt * swarm[i].v;
        boundary_condition(swarm[i].x);
    }
    force = cal_force(swarm);
    for (size_t i(0); i < N; ++i) {
        swarm[i].v = swarm[i].v + 0.5 * para.dt * swarm[i].mass_inv * force[i];
    }
}

size_t shuffle(vector<ptcl_t>& swarm) {
    // shuffle particles with Monte Carlo
    size_t Naccept(0);
    const size_t N(cal_N(swarm));
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
                dU += LJ(swarm[j], swarm[i]) - LJ(swarm[j], last_ptcl);
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
        ioer::info("c!");
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
            // adjust velocities for the rest particles
            vector<double> v_del(swarm[idx].v);
            vector<double> v_abs_sum(3, 0.0);
            for (int i(0); i < N; ++i) {
                if (i != idx) {
                    v_abs_sum = v_abs_sum + abs(swarm[i].v);
                }
            }
            for (int i(0); i < N; ++i) {
                if (i != idx) {
                    swarm[i].v = swarm[i].v + v_del * abs(swarm[i].v) / v_abs_sum;
                }
            }
            // remove the particle
            ioer::info("d!");
            swarm.erase(swarm.begin() + idx);
        }
    }
}

size_t MC_step(size_t istep, vector<ptcl_t>& swarm) {
    // single MC step
    size_t Naccept = shuffle(swarm);
    /*
    if (randomer::rand() < 0.5) {
        create(swarm);
    }
    else {
        destruct(swarm);
    }
    */
    return Naccept;
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
    ioer::info("# V = ", para.V, " Vc = ", para.Vc, " mu = ", para.mu);
    ioer::info("# Nstep_eql = ", para.Nstep_eql, " Nstep_run = ", para.Nstep_run);

    timer::tic();
    // init
    vector<ptcl_t> swarm(init_swarm(108));
    /*
    vector<ptcl_t> swarm {
        ptcl_t(
                vector<double>(3, 0.5 * para.L),
                randomer::maxwell_dist(para.mass, para.kT),
                para.mass,
                para.kT
              )
    };
    */

    // equilibrate
    for (size_t istep(0); istep < para.Nstep_eql; ++istep) {
        ioer::tabout(istep, cal_N(swarm) / para.V);
        MC_step(istep, swarm);
        //MD_step(istep, swarm);
    }

    // sampling
    double rho(0.0);
    double U_per_ptcl(0.0);
    double Utot(0.0);
    size_t Naccept(0);
    for (size_t istep(0); istep < para.Nstep_run; ++istep) {
        //rho += cal_N(swarm) / para.V;
        //U_per_ptcl += cal_U(swarm) / cal_N(swarm);
        size_t N(cal_N(swarm));
        Utot += cal_U(swarm) + LJ.U_LRC0 * N * N;
        //ioer::tabout(istep, rho / (istep + 1), U_per_ptcl / (istep + 1));
        ioer::tabout(istep, Utot / (istep + 1), Utot / (istep + 1) / N);

        Naccept += MC_step(istep, swarm);
        //MD_step(istep, swarm);
    }

    // output
    //ioer::tabout(rho / para.Nstep_run, U_per_ptcl / para.Nstep_run);
    ioer::info("# acc ratio = ", static_cast<double>(Naccept) / cal_N(swarm) / para.Nstep_run);
    ioer::info("# ", timer::toc());

    return 0;
}
