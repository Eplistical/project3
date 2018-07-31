#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <fstream>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "boost/timer/timer.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
using namespace std;

// config
struct Para {
    // basic
    double V = 512; // 555.55555555555555556;
    double L = pow(V, 1.0 / 3.0); 
    double rc = 2.5; //3;

    double kT = 1.0; // 0.9;
    double beta = 1.0 / kT;
    double mu = -3.2;

    // propagation
    uint32_t Nstep = 1e4;

    // MC related
    double dxmax = 0.1;
    bool dxmax_adaptive = true;
    bool prepinit = false;
    string conffile = "conf.dat";

    // others 
    uint32_t Anastep = 100;
    uint32_t random_seed = 0;

    public:
    void show() {
        ioer::keyval()
            ("# V", V)
            ("# L", L)
            ("# rc", rc)
            ("# kT", kT)
            ("# mu", mu)
            ("# Nstep", Nstep)
            ("# dxmax", dxmax)
            ("# dxmax_adaptive", dxmax_adaptive)
            ("# prepinit", prepinit)
            ("# conffile", conffile)
            ("# Anastep", Anastep)
            ("# random_seed", random_seed)
            ;
    }
} para;

// arg parser

bool argparse(int argc, char** argv, bool output_flag = true)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("kT", po::value<double>(&para.kT), "temperature")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("rc", po::value<double>(&para.rc), "cutoff range")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("dxmax_adaptive", po::value<bool>(&para.dxmax_adaptive), "if true, the program update dxmax according to acceptance ratio")
        ("Nstep", po::value<uint32_t>(&para.Nstep), "time step")
        ("Anastep", po::value<uint32_t>(&para.Anastep), "analysis step interval")
        ("conffile", po::value<string>(&para.conffile), "file for configuration")
        ("prepinit", po::value<bool>(&para.prepinit), "if true, prepare initial configuration")
        ("random_seed", po::value<uint32_t>(&para.random_seed), "random seed")
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
    if (vm.count("kT")) {
        para.beta = 1.0 / para.kT;
    }
    return true;
}

// functions

void read_conf(vector<double>& x, const string& conffile) {
    ioer::h5file_t in(conffile, ios::in);
    uint32_t N;
    in.read_attr("para", "N", N);
    in.read_dataset("x", x);
    assert(N * 3 == x.size());
    in.close();
}

void write_conf(const vector<double>& x, const string& conffile) {
    ioer::h5file_t out(conffile, ios::out);
    uint32_t N(x.size() / 3);
    out.create_dataset("para", vector<double>{0.0});
    out.create_attr("para", "N", N);
    out.create_dataset("x", x);
    out.close();
}

inline void raw_pair_energy(const double r2, double& U, double& W) {
    //assert(r2 > 0.5);
    double ir2, ir6;
    ir2 = 1.0 / r2;
    ir6 = ir2 * ir2 * ir2;
    U = 4.0 * ir6 * (ir6 - 1.0);
    W = 16.0 * ir6 * (ir6 - 0.5);
}

inline void pair_energy(const double* x1, const double* x2,
        const double rc, const double Urc, const double L,
        double& U, double& W)
{
    // caluclate pair energy between x1 & x2
    const double rc2(rc * rc);
    double dxk, r2;
    r2 = 0.0;
    for (int k(0); k < 3; ++k) {
        dxk = x1[k] - x2[k];
        // dxk -= L * round(dxk / L);

        if (dxk > 0.5 * para.L) dxk -= para.L;
        else if (dxk < -0.5 * para.L) dxk += para.L;

        assert(dxk >= -0.5 * para.L and dxk <= 0.5 * para.L);
        r2 += dxk * dxk;
    }

    if (r2 < rc2) {
        raw_pair_energy(r2, U, W);
        U -= Urc;
    }
    else {
        U = 0.0;
        W = 0.0;
    }

    /*
    double sr2, sr6;
    if (r2 < rc2) {
        sr2 = 1.0 / r2;
        sr6 = sr2 * sr2 * sr2;
        U = 4.0 * sr6 * (sr6 - 1.0);
        W = 16.0 * sr6 * (sr6 - 0.5);
    }
    else {
        U = 0.0;
        W = 0.0;
    }
    */
    /*
        U = 0.0;
        W = 0.0;
    */
}


void one_energy(const vector<double>& x, const uint32_t ofs, 
        const double rc, const double Urc, const double L,
        double& U, double& W) 
{
    // calculate the U & W of given atom w/ all other atoms
    if (not x.empty()) {
        const uint32_t _3N(x.size());

        U = 0.0;
        W = 0.0;

        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
            if (ofs_i != ofs) {
                pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij);
                U += Uij;
                W += Wij;
            }
        }
    }
    else {
        U = 0.0;
        W = 0.0;
    }
}


void all_energy(const vector<double>& x, 
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0,
        double& U, double& W) 
{
    // calculate the sum of U & W of all energies
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        const uint32_t _3N_minus_3(x.size() - 3);
        const uint32_t N(_3N / 3);

        U = ULRC0 * N * N;
        W = WLRC0 * N * N;

        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N_minus_3; ofs_i += 3) {
            for (uint32_t ofs_j(ofs_i + 3); ofs_j < _3N; ofs_j += 3) {
                pair_energy(&x[ofs_i], &x[ofs_j], rc, Urc, L, Uij, Wij);
                U += Uij;
                W += Wij;
            }
        }
    }
    else {
        U = 0.0;
        W = 0.0;
    }
}

void potin(const vector<double>& newx, const vector<double>& x, 
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for a new particle newx
    const uint32_t _3N(x.size());
    const uint32_t N(_3N / 3);

    dU = (2.0 * static_cast<double>(N) + 1.0) * ULRC0;
    dW = (2.0 * static_cast<double>(N) + 1.0) * WLRC0;

    double Uij, Wij;
    for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
        pair_energy(&newx[0], &x[ofs_i], rc, Urc, L, Uij, Wij);
        dU += Uij;
        dW += Wij;
    }
}

void potout(const vector<double>& x, const uint32_t ofs,
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for removing an existing particle x[ofs]
    const uint32_t _3N(x.size());
    const uint32_t N(_3N / 3);
    assert(N > 0);

    dU = -(2.0 * static_cast<double>(N) - 1.0) * ULRC0;
    dW = -(2.0 * static_cast<double>(N) - 1.0) * WLRC0;

    double Uij, Wij;
    for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
        if (ofs_i != ofs) {
            pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij);
            dU -= Uij;
            dW -= Wij;
        }
    }
}

inline bool decide(double x) {
    // generate a prob by e^-x
    if (x < 75.0) {
        if (x < 0.0) {
            return true;
        }
        else if (randomer::rand() < exp(-x)) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}

void muVTMC() 
{
    // shifted energy
    double Urc(0.0), Wrc(0.0);
    raw_pair_energy(para.rc * para.rc, Urc, Wrc);

    // long range correction
    double sr3, sr9, ULRC0(0.0), WLRC0(0.0);
    sr3 = pow(1.0 / para.rc, 3);
    sr9 = sr3 * sr3 * sr3;
    /*
    ULRC0 = 8.0 / 9.0 * M_PI / para.V * (sr9 - 3.0 * sr3);
    WLRC0 = 32.0 / 9.0 * M_PI / para.V * (sr9 - 1.5 * sr3);
    */

    // prepare/load init configuration
    vector<double> x;
    if (para.prepinit) {
        x = randomer::vrand(3 * 1, -0.5 * para.L, 0.5 * para.L);
    }
    else {
        read_conf(x, para.conffile);
    }
    ioer::info("# init N = ", x.size() / 3);

    // calc init U & W
    double U, W;
    all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, U, W);
    ioer::info("# init U = ", U, " init W = ", W);

    // main MC part
    vector<double> oldx(3), newx(3);
    uint32_t _3N, N;
    double Uold, Wold, Unew, Wnew;
    double dU, dW, beta_dU, dCB, dDB;
    double Usum(0.0), Psum(0.0), rhosum(0.0); 
    uint32_t Naccept(0), Nmove(0);
    uint32_t Nsamp(0);
    double acc_ratio;

    // for debug
    double Utot, Wtot;

    // loop over steps
    ioer::info("# start looping ... ");
    for (uint32_t istep(0); istep < para.Nstep; ++istep) {

        // shuffle, loop over atoms
        _3N = x.size();
        N = _3N / 3;

        for (uint32_t ofs(0); ofs < _3N; ofs += 3) {

            one_energy(x, ofs, para.rc, Urc, para.L, Uold, Wold);
            for (int k(0); k < 3; ++k) {
                oldx[k] = x[ofs + k];
                x[ofs + k] += randomer::rand(-para.dxmax, para.dxmax);
                x[ofs + k] -= para.L * round(x[ofs + k] / para.L);
            }
            one_energy(x, ofs, para.rc, Urc, para.L, Unew, Wnew);

            dU = Unew - Uold;
            dW = Wnew - Wold;
            beta_dU = para.beta * dU;

            //if (decide(beta_dU)) {
            if (randomer::rand() < exp(-beta_dU)) {
                U += dU;
                W += dW;
                Naccept += 1;
                //ioer::info("s: dU = ", dU);
            }
            else {
                copy(oldx.begin(), oldx.end(), x.begin() + ofs);
            }
            Nmove += 1;

        } // shuffle, loop over atoms


        // create or destruct
        if (randomer::rand() < 0.5) {
            // create
            newx = randomer::vrand(3, -0.5 * para.L, 0.5 * para.L);
            potin(newx, x, para.rc, Urc, para.L, ULRC0, WLRC0, dU, dW);

            /*
            dCB = para.beta * dU - para.beta * para.mu - log(para.V / (N + 1));
            if (decide(dCB)) {
            */
            const double prob(exp(-para.beta * dU + para.beta * para.mu) * para.V / (N + 1));
            if (randomer::rand() < prob) {
                //ioer::info("create: N = ", N, " dU = ", dU, " prob = ", prob, " newx = ", newx);
                N += 1;

                x.resize(_3N + 3);
                copy(newx.begin(), newx.end(), x.end() - 3);
                U += dU;
                W += dW;
            }
        }
        else {
            // destruct
            if (N > 0) {
                uint32_t ofs(randomer::choice(N) * 3);
                potout(x, ofs, para.rc, Urc, para.L, ULRC0, WLRC0, dU, dW);

                /*
                dDB = para.beta * dU + para.beta * para.mu - log(N / para.V);
                if (decide(dDB)) {
                */
                const double prob(exp(-para.beta * dU - para.beta * para.mu) * N / para.V);
                if (randomer::rand() < prob) {
                    //ioer::info("destruct: N = ", N, " dU = ", dU, " prob = ", prob);
                    N -= 1;

                    copy(x.end() - 3, x.end(), x.begin() + ofs);
                    x.resize(_3N - 3);
                    U += dU;
                    W += dW;

                    //double Utot, Wtot;
                    //all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, Utot, Wtot);
                    //ioer::info("tot U = ", Utot, " U = ", U);
                    //ioer::info("d: dU = ", dU);
                }
            }
        }

        assert(N < 1000);
        for (double xi : x) {
            assert(xi <= 0.5 * para.L and xi >= -0.5 * para.L);
        }
        
        // statistics
        Nsamp += 1;
        if (N > 0) {
            rhosum += N / para.V;
            Usum += U / N;
            Psum += N / para.V * para.kT + W  / para.V;
        }

        // Analyze
        if (istep % para.Anastep == 0) {

            //double Utot, Wtot; all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, Utot, Wtot);

            // print statistics
            ioer::tabout(Nsamp, x.size() / 3, rhosum / Nsamp, Usum / Nsamp, Psum / Nsamp);
            // save conf
            write_conf(x, para.conffile);
            // adjust dxmax
            if (para.dxmax_adaptive) {
                acc_ratio = static_cast<double>(Naccept) / Nmove;
                if (acc_ratio < 0.4) {
                    para.dxmax *= 0.95;
                }
                else if (acc_ratio > 0.6) {
                    para.dxmax *= 1.05;
                }
                Nmove = 0;
                Naccept = 0;
            }
        }

    } // loop over steps

        /*
    N = x.size() / 3;
    double uij, wij, utot(0.0);
    vector<double> dx;
    for (int i(0); i < N; ++i) {
        vector<double> xi(x.begin() + i*3, x.begin() + i*3 + 3);
        for (int j(0); j < N; ++j) {
            if (i < j) {
                vector<double> xj(x.begin() + j*3, x.begin() + j*3 + 3);

                pair_energy(&xi[0], &xj[0], para.rc, Urc, para.L, uij, wij);

                utot += uij;

                dx = xi - xj;
                for (double& dxk : dx) {
                    if (dxk > 0.5 * para.L) 
                        dxk -= para.L;
                    else if (dxk < -0.5 * para.L)
                        dxk += para.L;
                }

                ioer::tabout(i, j, 
                        xi - xj,
                        dx,
                        norm(dx),
                        uij, utot
                        );
            }
        }
    }
    */

} // main


int main(int argc, char** argv) {
    /*
    // test
    vector<double> x1(3, 0.0), x2(3, 0.0), x3(3, 0.0);
    double Uij, Wij;

    for (double xi(0.7); xi < para.L - 0.7; xi += 0.001) {
        x2[0] = xi;
        pair_energy(&x1[0], &x2[0], para.rc, Urc, para.L, Uij, Wij);
        ioer::tabout(xi, Uij, Wij);
    }

    return 0;
    */

    if (argparse(argc, argv) == false) {
        return 0;
    }
    else {
        para.show();
        randomer::seed(para.random_seed);

        boost::timer::cpu_timer cpu_timer;
        muVTMC();
        ioer::info("# ", cpu_timer.format(4));
    }
    return 0;
}
