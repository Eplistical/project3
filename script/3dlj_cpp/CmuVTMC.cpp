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
    double V = 512;
    double L = pow(V, 1.0 / 3.0); 
    double rc = 2.5;

    double kT = 1.0; 
    double mu = -3.2;

    // propagation
    uint32_t Nstep = 1e4;

    // MC related
    double dxmax = 0.1;
    bool prepinit = false;
    uint32_t N0 = 108;
    string conffile = "conf.dat";

    // others 
    uint32_t Anastep = 100;
    uint32_t random_seed = 0;

    public:
    void show(ioer::output_t& out) {
        out.info("# --- CONFIG PARAMETERS --- ");
        out.keyval()
            ("# V", V)
            ("# L", L)
            ("# rc", rc)
            ("# kT", kT)
            ("# mu", mu)
            ("# Nstep", Nstep)
            ("# dxmax", dxmax)
            ("# prepinit", prepinit)
            ("# N0", N0)
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
        ("V", po::value<double>(&para.V), "box volume")
        ("kT", po::value<double>(&para.kT), "temperature")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("rc", po::value<double>(&para.rc), "cutoff range")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("Nstep", po::value<uint32_t>(&para.Nstep), "time step")
        ("Anastep", po::value<uint32_t>(&para.Anastep), "analysis step interval")
        ("conffile", po::value<string>(&para.conffile), "file for configuration")
        ("prepinit", po::value<bool>(&para.prepinit), "if true, prepare initial configuration")
        ("N0", po::value<uint32_t>(&para.N0), "initial # of atoms to prepare")
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
    if (vm.count("V")) {
        para.L = pow(para.V, 1.0 / 3.0);
    }
    return true;
}

// configuration I/O functions

void read_conf(vector<double>& x, const string& conffile) 
{
    ioer::h5file_t in(conffile, ios::in);
    uint32_t N;
    in.read_attr("para", "N", N);
    in.read_dataset("x", x);
    assert(N * 3 == x.size());
    in.close();
}

inline void write_conf(const vector<double>& x, const string& conffile) 
{
    ioer::h5file_t out(conffile, ios::out);
    uint32_t N(x.size() / 3);
    out.create_dataset("para", vector<double>{0.0});
    out.create_attr("para", "N", N);
    out.create_dataset("x", x);
    out.close();
}

// energy functions

inline void raw_pair_energy(const double r2, double& U, double& W) 
{
    // calc pair energy for r2
    double ir2, ir6;
    ir2 = 1.0 / r2;
    ir6 = ir2 * ir2 * ir2;
    U = 4.0 * ir6 * (ir6 - 1.0);
    W = 16.0 * ir6 * (ir6 - 0.5);
}

inline void pair_energy(const double* x1, const double* x2,
        const double rc, const double Urc, const double L,
        double& U, double& W, bool& invin)
{
    // calc pair energy between x1 & x2
    // subject to min imag & cut-off
    const double rc2(rc * rc);
    double dxk, r2;
    r2 = 0.0;
    for (int k(0); k < 3; ++k) {
        dxk = x1[k] - x2[k];
        dxk -= L * round(dxk / L);
        r2 += dxk * dxk;
    }

    if (r2 < rc2) {
        invin = true;
        raw_pair_energy(r2, U, W);
        U -= Urc;
    }
    else {
        invin = false;
        U = 0.0;
        W = 0.0;
    }
}

inline void one_energy(const vector<double>& x, const uint32_t ofs, 
        const double rc, const double Urc, const double L,
        double& U, double& W, uint32_t& Nvin) 
{
    // calc U & W for a given atom w/ other atoms
    bool invin;
    Nvin = 0;
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
            if (ofs_i != ofs) {
                pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij, invin);
                if (invin) {
                    Nvin += 1;
                    U += Uij;
                    W += Wij;
                }
            }
        }
    }
}

void all_energy(vector<double>& x, 
        const double rc, const double Urc, const double L, 
        const double ULRC0, const double WLRC0, 
        double& U, double& W, uint32_t& Nvin)
{
    // calc total U & W for a given x
    bool invin;
    Nvin = 0;
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        const uint32_t _3N_3(_3N - 3);
        const uint32_t N(_3N / 3);
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N_3; ofs_i += 3) {
            for (uint32_t ofs_j(ofs_i + 3); ofs_j < _3N; ofs_j += 3) {
                pair_energy(&x[ofs_i], &x[ofs_j], rc, Urc, L, Uij, Wij, invin);
                if (invin) {
                    Nvin += 1;
                    U += Uij;
                    W += Wij;
                }
            }
        }
        U += ULRC0 * N * N;
        W += WLRC0 * N * N;
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
    uint32_t Nvin;
    bool invin;

    Nvin = 0;
    dU = 0.0;
    dW = 0.0;

    double Uij, Wij;
    for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
        pair_energy(&newx[0], &x[ofs_i], rc, Urc, L, Uij, Wij, invin);
        if (invin) {
            Nvin += 1;
            dU += Uij;
            dW += Wij;
        }
    }

    dU += (2.0 * static_cast<double>(N) + 1.0) * ULRC0;
    dW += (2.0 * static_cast<double>(N) + 1.0) * WLRC0;
}

void potout(const vector<double>& x, const uint32_t ofs,
        const double rc, const double Urc, const double L,
        const double ULRC0, const double WLRC0, 
        double& dU, double& dW)
{
    // calc dU & dW for removing an existing particle x[ofs]
    assert(not x.empty());
    const uint32_t N(x.size() / 3);
    uint32_t Nvin;

    Nvin = 0;
    dU = 0.0;
    dW = 0.0;

    one_energy(x, ofs, rc, Urc, L, dU, dW, Nvin);

    dU = -dU - (2.0 * static_cast<double>(N) - 1.0) * ULRC0;
    dW = -dW - (2.0 * static_cast<double>(N) - 1.0) * WLRC0;
}

void tail_correction(const double rc, const double V,
        double& Urc, double& ULRC0, double& WLRC0) 
{
    // calcualte tail correction energy
    ULRC0 = 0.0;
    WLRC0 = 0.0;
    Urc = 0.0;

    const bool lrc_flag(false);
    const bool shift_flag(true); // truncated + shifted potential

    if (lrc_flag) {
        // long range correction
        double irc3, irc9;
        irc3 = pow(1.0 / rc, 3);
        irc9 = irc3 * irc3 * irc3;
        ULRC0 = 8.0 / 9.0 * M_PI / V * (irc9 - 3.0 * irc3);
        WLRC0 = 32.0 / 9.0 * M_PI / V * (irc9 - 1.5 * irc3);
    }

    if (shift_flag) {
        // shifted energy
        double Wrc, rc2(rc * rc);
        raw_pair_energy(rc2, Urc, Wrc);
    }
}

// configuration move functions

inline bool decide(double x) {
    // generate a prob by e^-x
    if (x > 75.0) 
        return false;
    else if (x < 0.0)
        return true;
    else if (randomer::rand() < exp(-x))
        return true;
    else 
        return false;
}

bool shuffle(vector<double>& x, const uint32_t ofs,
        double& U, double& W,
        const double kT, const double dxmax, const double L,
        const double rc, const double Urc)
{
    // shuffle the particle x[ofs:ofs+3] w/ MC algorithm
    uint32_t _3N = x.size();
    assert (ofs % 3 == 0 and ofs < _3N);

    uint32_t Nvinold, Nvinnew;
    double Uold, Wold, Unew, Wnew;
    double dU, dW;
    vector<double> oldx(3);
    one_energy(x, ofs, rc, Urc, L, Uold, Wold, Nvinold);

    for (int k(0); k < 3; ++k) {
        oldx[k] = x[ofs + k];
        x[ofs + k] += randomer::rand(-dxmax, dxmax);
        x[ofs + k] -= L * round(x[ofs + k] / L); // periodic condition
    }
    
    one_energy(x, ofs, rc, Urc, L, Unew, Wnew, Nvinnew);

    dU = Unew - Uold;
    dW = Wnew - Wold;

    if (decide(dU / kT)) {
        U += dU;
        W += dW;
        return true;
    }
    else {
        copy(oldx.begin(), oldx.end(), x.begin() + ofs);
        return false;
    }
}

bool create(vector<double>& x, const vector<double>& newx,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    uint32_t N(_3N / 3);
    double dU, dW;
    double dCB;
    potin(newx, x, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dCB = dU / kT - mu / kT  - log(V / (N + 1));
    if (decide(dCB)) {
        x.resize(_3N + 3);
        copy(newx.begin(), newx.begin() + 3, x.end() - 3);
        U += dU;
        W += dW;
        return true;
    }
    return false;
}

bool destruct(vector<double>& x, const uint32_t ofs,
        double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    assert(ofs % 3 == 0 and ofs < _3N);

    uint32_t N(_3N / 3);
    double dU, dW;
    double dDB;
    potout(x, ofs, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dDB = dU / kT + mu / kT - log(N / V);
    if (decide(dDB)) {
        copy(x.end() - 3, x.end(), x.begin() + ofs);
        x.resize(_3N - 3);
        U += dU;
        W += dW;
        return true;
    }
    return false;
}

// main functions

void muVTMC() 
{
    // output config
    ioer::output_t out("");
    out.set_precision(4);

    // show para
    para.show(out);
    out.info("# --- PROGRAM BEGINS --- ");
    out.info("# exe name: muVTMC ");

    // tail corrections
    double ULRC0, WLRC0, Urc;
    tail_correction(para.rc, para.V, Urc, ULRC0, WLRC0);
    out.info("# tail correction: Urc = ", Urc, " ULRC0 = ", ULRC0, " WLRC0 = ", WLRC0);

    // init configuration
    vector<double> x;
    double U, W;
    uint32_t Nvin;
    if (para.prepinit) {
        x = randomer::vrand(3 * para.N0, -0.5 * para.L, 0.5 * para.L);
    }
    else {
        read_conf(x, para.conffile);
    }
    all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, U, W, Nvin);
    out.info("# init configuration: N = ", x.size() / 3, " init U = ", U, " init W = ", W);

    // main MC part
    double randnum;
    const double shuffle_frac(0.75), create_frac(1.0 - (1.0 - shuffle_frac) * 0.5);
    uint32_t N;
    uint32_t Naccept(0), Nmove(0);
    double Usum(0.0), Psum(0.0), rhosum(0.0);
    uint32_t Nsamp(0);
    double Unow, Pnow, rhonow;
    out.info("# start MC looping ... ");
    out.info("# shuffle frac = ", shuffle_frac, " create frac = ", create_frac);
    out.tabout("# Nsamp", "rho", "U/N", "P", "<rho>", "<U/N>", "<P>");
    for (uint32_t istep(0); istep < para.Nstep; ++istep) {
        N = x.size() / 3;

        randnum = randomer::rand();
        if (randnum < shuffle_frac) {
            uint32_t ofs(randomer::choice(N) * 3);
            Naccept += shuffle(x, ofs, U, W, para.kT, para.dxmax, para.L, para.rc, Urc);
        }
        else if (randnum < create_frac) {
            vector<double> newx(randomer::vrand(3, -0.5 * para.L, 0.5 * para.L));
            Naccept += create(x, newx, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
        }
        else {
            uint32_t ofs(randomer::choice(N) * 3);
            Naccept += destruct(x, ofs, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
        }
        Nmove += 1;

        // statistics
        Nsamp += 1;
        N = x.size() / 3;
        if (N > 0) {
            rhonow = N / para.V;
            Unow = U / N;
            Pnow = N / para.V * para.kT + W / para.V;

            rhosum += rhonow;
            Usum += Unow;
            Psum += Pnow;
        }
        // output & save
        if (istep % para.Anastep == 0) {
            // print statistics & save configure
            out.tabout(Nsamp, 
                    rhonow, Unow, Pnow,
                    rhosum / Nsamp,
                    Usum / Nsamp,
                    Psum / Nsamp
                    );
            write_conf(x, para.conffile);
        }
    } 
    out.tabout(Nsamp, 
            rhonow, Unow, Pnow,
            rhosum / Nsamp,
            Usum / Nsamp,
            Psum / Nsamp
            );
    write_conf(x, para.conffile);

    //final output
    out.drawline('\n', 1);
    out.tabout("# <rho> = ", rhosum / Nsamp);
    out.tabout("# <U/N> = ", Usum / Nsamp);
    out.tabout("# <p> = ", Psum / Nsamp);
    out.tabout("# acc ratio = ", static_cast<double>(Naccept) / Nmove);
    out.drawline('\n', 1);

    out.info("# --- PROGRAM END --- ");
} 

// main

int main(int argc, char** argv) {

    if (argparse(argc, argv) == false) {
        return 0;
    }
    else {
        randomer::seed(para.random_seed);

        boost::timer::cpu_timer cpu_timer;
        muVTMC();
        ioer::info("# ", cpu_timer.format(4));
    }
    return 0;
}
