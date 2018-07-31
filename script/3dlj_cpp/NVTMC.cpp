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
 * NVT MC simulation for 3D LJ fluid
 */

// config
struct Para {
    // basic
    const size_t N = 108;
    const double V = N / 0.9;
    const double L = pow(V, 1.0 / 3.0); 
    const double rc = 3.0;

    double kT = 0.9;
    double beta = 1.0 / kT;
    double mass = 1.0;
    size_t random_seed = 1287234;

    // propagation
    size_t Nstep_eql = 1e5;
    size_t Nstep_run = 7e5;

    // MC related
    double dxmax = 0.5;
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
            // W_LRC0 = W_LRC / N**2
            W_LRC0 = 32.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 1.5 * pow(rc, -3));
            // mu_LRC0 = mu_LRC / N
            mu_LRC0 = 16.0 / 9.0 * M_PI / Vbox * (pow(rc, -9) - 3.0 * pow(rc, -3));
        }

    public:
        double operator()(const double* x1, const double* x2) {
            /**
             * calculate shifted potential, virial between coordinates 1 & 2
             * return shifted potential
             */
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

            // calc U_shifted, virial
            const double r6(r2 * r2 * r2);
            const double r6_inv(1.0 / r6);
            U_shifted = 4.0 * r6_inv * (r6_inv - 1.0) - U_rc;
            virial = 48.0 * r6_inv * (r6_inv - 0.5);

            return U_shifted;
        }

    public:
        double U_shifted;
        double virial;

    public:
        double rc, rc2, U_rc;
        double Lbox, Vbox;
        double U_LRC0, W_LRC0, mu_LRC0;
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
void cal_U_and_W(const vector<double>& x, double& U, double& W) {
    const size_t _3N(para.N * 3);
    const size_t _3N_minus_3(_3N - 3);
    U = 0.0;
    W = 0.0;

    for (size_t offset_i(0); offset_i < _3N_minus_3; offset_i += 3) {
        for (size_t offset_j(offset_i + 3); offset_j < _3N; offset_j += 3) {
            LJ(&x[offset_i], &x[offset_j]);
            U += LJ.U_shifted;
            W += LJ.virial;
        }
    }

    W /= 3.0;

    // long range correction
    U += LJ.U_LRC0 * para.N * para.N;
    W += LJ.W_LRC0 * para.N * para.N;
}

size_t shuffle(vector<double>& x) {
    // shuffle particles with Monte Carlo
    const size_t _3N(3 * para.N);

    size_t Naccept(0);
    vector<double> dx(randomer::vrand(para.N * 3, -para.dxmax, para.dxmax));
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
        const double beta_dU(para.beta * dU);
        bool acc_flag(false);
        if (beta_dU >= 100.0) {
            acc_flag = false;
        }
        else if (beta_dU <= 0.0) {
            acc_flag = true;
        }
        else {
            acc_flag = (randomer::rand() < exp(-beta_dU));
        }
        // move
        if (acc_flag) {
            copy(new_x.begin(), new_x.end(), x.begin() + offset_i);
            Naccept += 1;
        }
    }
    return Naccept;
}

int main(int argc, char** argv) {
    // parse args
    if (argparse(argc, argv) == false) {
        return 0;
    }
    timer::tic();

    // setup
    randomer::seed(para.random_seed);
    vector<double> x(randomer::vrand(para.N * 3, 0.0, para.L));

    // equilibrate
    for (size_t istep(1); istep <= para.Nstep_eql; ++istep) {
        shuffle(x);
    }

    // sampling
    double Usum(0.0);
    double Wsum(0.0);
    double W, U;
    size_t Naccept(0);
    for (size_t istep(1); istep <= para.Nstep_run; ++istep) {
        cal_U_and_W(x, U, W);
        Usum += U;
        Wsum += W;

        Naccept += shuffle(x);
    }

    // post-processing
    double avgU(Usum / para.Nstep_run / para.N);

    double avgP((para.N * para.kT + Wsum / para.Nstep_run) / para.V);

    double acc_ratio(static_cast<double>(Naccept) / para.N / para.Nstep_run);

    // output
    ioer::output_t out("");
    out.set_precision(8);
    out.info("# rc = ", para.rc, " kT = ", para.kT);
    out.info("# V = ", para.V, " rho = ", para.N / para.V, " N = ", para.N);
    out.info("# Nstep_eql = ", para.Nstep_eql, " Nstep_run = ", para.Nstep_run);
    out.info("# random_seed = ", para.random_seed);

    out.info("# avg U = ", avgU);
    out.info("# avg P = ", avgP);
    out.info("# acc ratio = ", acc_ratio);

    out.info("# ", timer::toc());

    return 0;
}
