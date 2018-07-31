#include <cmath>
#include <cassert>
#include <algorithm>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "misc/timer.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
using namespace std;

/**
 * muVT MC simulation for 3D LJ fluid
 */

// config
struct Para {
    // basic
    const double V = 512;
    const double L = pow(V, 1.0 / 3.0); 
    const double rc = 2.5;
    const double r_overlap = 0.5;

    double mu = -3.2;
    double kT = 1.0;
    double beta = 1.0 / kT;
    double mass = 1.0;
    size_t random_seed = 1287234;

    // propagation
    size_t Nstep_eql = 1e5;
    size_t Nstep_run = 4e5;

    // MC related
    double dxmax = 0.5;
} para;

// potential
struct LJ_shifted {
    public:
        LJ_shifted(double RC, double LBOX, double R_OVERLAP) noexcept 
            : rc(RC), Lbox(LBOX), r_overlap(R_OVERLAP)
        {
            // -- basic parameters -- //
            rc2 = rc * rc;
            r_overlap2 = r_overlap * r_overlap;
            Vbox = Lbox * Lbox * Lbox;

            // -- shifted energy -- //
            U_rc = 4.0 * (pow(rc, -12) - pow(rc, -6));

            // -- long range correction -- //
            
            // U_LRC0 = U_LRC / N**2, W_LRC0 = W_LRC / N**2, mu_LRC0 = mu_LRC / N
            U_LRC0 = 8.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 3.0 * pow(rc, -3)); 
            W_LRC0 = 32.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 1.5 * pow(rc, -3));
            mu_LRC0 = 16.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 3.0 * pow(rc, -3));
            /*
            U_LRC0 = 0.0;
            W_LRC0 = 0.0;
            mu_LRC0 = 0.0;
            */
        }

    public:
        double operator()(const double* x1, const double* x2) {
            /**
             * calculate shifted potential, virial between coordinates 1 & 2
             * return shifted potential
             */
            overlap = false;
            double r2(0.0);
            double dxi;
            for (int i(0); i < 3; ++i) {
                dxi = x1[i] - x2[i];
                // min imag
                dxi -= Lbox * round(dxi / Lbox);
                assert(dxi <= 0.5 * Lbox and dxi >= -0.5 * Lbox);

                // get r2
                if (abs(dxi) > rc) {
                    U_shifted = 0.0;
                    virial = 0.0;
                    return U_shifted;
                }
                else {
                    r2 += dxi * dxi;
                    if (r2 > rc2) {
                        U_shifted = 0.0;
                        virial = 0.0;
                        return U_shifted;
                    }
                }
            }
            // determine overlap
            if (r2 < r_overlap2) {
                overlap = true;
            }

            // calc U_shifted, virial
            const double r6(r2 * r2 * r2);
            const double r6_inv(1.0 / r6);
            U_shifted = 4.0 * r6_inv * (r6_inv - 1.0) - U_rc;
            virial = -48.0 * r6_inv * (r6_inv - 0.5);

            return U_shifted;
        }

    public:
        double U_shifted;
        double virial;

    public:
        double rc, rc2, U_rc;
        double r_overlap, r_overlap2;
        double Lbox, Vbox;
        double U_LRC0, W_LRC0, mu_LRC0;
        bool overlap;
} LJ(para.rc, para.L, para.r_overlap);

// arg parser

bool argparse(int argc, char** argv, bool output_flag = true)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("Nstep_eql", po::value<size_t>(&para.Nstep_eql), "time step for equilibrating the system")
        ("Nstep_run", po::value<size_t>(&para.Nstep_run), "time step for collecting data")
        ("random_seed", po::value<size_t>(&para.random_seed), "random seed")
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

// functions

size_t shuffle(vector<double>& x) {
    // shuffle particles with Monte Carlo
    const size_t _3N(x.size());

    size_t Naccept(0);
    vector<double> dx(randomer::vrand(_3N, -para.dxmax, para.dxmax));
    vector<double> new_x(3, 0.0);

    for (size_t offset_i(0); offset_i < _3N; offset_i += 3) {
        // move w/ boundary condition
        for (int k(0); k < 3; ++k) {
            new_x[k] = x[offset_i + k] + dx[offset_i + k];

            if (new_x[k] < 0.0) {
                new_x[k] += para.L;
            }
            else if (new_x[k] > para.L) {
                new_x[k] -= para.L;
            }
            assert (new_x[k] >= 0.0 and new_x[k] <= para.L);
        }
        // calc U
        double dU(0.0);
        for (size_t offset_j(0); offset_j < _3N; offset_j += 3) {
            if (offset_i != offset_j) {
                dU += LJ(&new_x[0], &x[offset_j]) - LJ(&x[offset_i], &x[offset_j]);
            }
        }
        // make decision
        const double prob(exp(-para.beta * dU));
        if (randomer::rand() < prob) {
            copy(new_x.begin(), new_x.end(), x.begin() + offset_i);
            Naccept += 1;
        }
    }
    return Naccept;
}


bool create(vector<double>& x) {
    // create a new particle
    const size_t _3N(x.size());
    const size_t N(_3N / 3);
    // generate a new particle
    vector<double> new_x(randomer::vrand(3, 0.0, para.L));
    // calc dU
    double dU(0.0);
    for (size_t offset_i(0); offset_i < _3N; offset_i += 3) {
        dU += LJ(&new_x[0], &x[offset_i]);
    }
    dU += LJ.U_LRC0 * (2 * N + 1);
    // make decision
    const double prob(exp(-para.beta * (dU - para.mu)) * para.V / (N + 1));
    if (randomer::rand() < prob) {
        x.resize(_3N + 3);
        copy(new_x.begin(), new_x.end(), x.end() - 3);
        return true;
    }
    else {
        return false;
    }
}


bool destruct(vector<double>& x) {
    // remove an existing particle
    const size_t _3N(x.size());
    const size_t N(_3N / 3);
    // if empty
    if (N == 0) {
        return false;
    }
    // pick the particle 
    size_t offset(randomer::choice(N) * 3);
    // calc dU
    double dU(0.0);
    for (size_t offset_i(0); offset_i < _3N; offset_i += 3) {
        if (offset_i != offset) {
            dU -= LJ(&x[offset], &x[offset_i]);
        }
    }
    dU -= LJ.U_LRC0 * (2 * N - 1);
    // make decision
    const double prob(exp(-para.beta * (dU + para.mu)) * N / para.V);
    if (randomer::rand() < prob) {
        copy(x.end() - 3, x.end(), x.begin() + offset);
        x.resize(_3N - 3);
        return true;
    }
    else {
        return false;
    }
}


void cal_U_and_W(const vector<double>& x, double& U, double& W) {
    const size_t _3N(x.size());
    const size_t _3N_minus_3(x.size() - 3);
    const size_t N(_3N / 3);
    const double rho(N / para.V);
    U = 0.0;
    W = 0.0;

    for (size_t offset_i(0); offset_i < _3N_minus_3; offset_i += 3) {
        for (size_t offset_j(offset_i + 3); offset_j < _3N; offset_j += 3) {
            LJ(&x[offset_i], &x[offset_j]);
            U += LJ.U_shifted;
            W += LJ.virial;
        }
    }

    W /= -3.0;
    
    // long range correction
    U += LJ.U_LRC0 * N * N;
    W += LJ.W_LRC0 * N * N;
}


int main(int argc, char** argv) {
    // parse args
    if (argparse(argc, argv) == false) {
        return 0;
    }
    timer::tic();

    // setup
    randomer::seed(para.random_seed);
    vector<double> x(randomer::vrand(3 * 1, 0.0, para.L));

    // equilibrate
    for (size_t istep(1); istep <= para.Nstep_eql; ++istep) {
        shuffle(x);
        if (randomer::rand() < 0.5) {
            create(x);
        }
        else {
            destruct(x);
        }
        ioer::tabout(istep, x.size() / 3.0 / para.V);
    }

    // sampling
    double Usum(0.0);
    double W, U;
    double Psum(0.0);
    double Nsum(0.0);
    double rhosum(0.0);

    size_t Naccept(0);
    size_t Nattempt(0);
    size_t Naccept_c(0);
    size_t Naccept_d(0);
    for (size_t istep(1); istep <= para.Nstep_run; ++istep) {
        size_t N(x.size()  / 3);
        Nsum += N;
        rhosum += N / para.V;

        cal_U_and_W(x, U, W);
        Usum += U / N;
        Psum += N * para.kT / para.V + W / para.V;

        Naccept += shuffle(x);
        if (randomer::rand() < 0.5) {
            Naccept_c += create(x);
        }
        else {
            Naccept_d += destruct(x);
        }
        
        ioer::tabout(istep, 
                rhosum / istep, 
                Usum / istep,
                Psum / istep,
                static_cast<double>(Naccept) / Nsum,
                LJ.U_LRC0,
                LJ.W_LRC0
                );
    }

    // post-processing

    // output
    ioer::output_t out("");
    out.set_precision(8);
    out.info("# rc = ", para.rc, " kT = ", para.kT, " mu = ", para.mu);
    out.info("# V = ", para.V, " L = ", para.L);
    out.info("# Nstep_eql = ", para.Nstep_eql, " Nstep_run = ", para.Nstep_run);
    out.info("# random_seed = ", para.random_seed);

    /*
    out.info("# avg U = ", avgU);
    out.info("# avg P = ", avgP);
    out.info("# acc ratio = ", acc_ratio);
    */

    out.info("# ", timer::toc());

    return 0;
}
