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
    void show() {
        ioer::keyval()
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

// functions

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
        double& U, double& W)
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
        raw_pair_energy(r2, U, W);
        U -= Urc;
    }
    else {
        U = 0.0;
        W = 0.0;
    }
}

inline void one_energy(const vector<double>& x, const uint32_t ofs, 
        const double rc, const double Urc, const double L,
        double& U, double& W) 
{
    // calc U & W for a given atom w/ other atoms
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N; ofs_i += 3) {
            if (ofs_i != ofs) {
                pair_energy(&x[ofs], &x[ofs_i], rc, Urc, L, Uij, Wij);
                U += Uij;
                W += Wij;
            }
        }
    }
}

void all_energy(vector<double>& x, 
        const double rc, const double Urc, const double L, 
        const double ULRC0, const double WLRC0, 
        double& U, double& W)
{
    // calc total U & W for a given x
    U = 0.0;
    W = 0.0;
    if (not x.empty()) {
        const uint32_t _3N(x.size());
        const uint32_t _3N_3(_3N - 3);
        const uint32_t N(_3N / 3);
        double Uij, Wij;
        for (uint32_t ofs_i(0); ofs_i < _3N_3; ofs_i += 3) {
            for (uint32_t ofs_j(ofs_i + 3); ofs_j < _3N; ofs_j += 3) {
                pair_energy(&x[ofs_i], &x[ofs_j], rc, Urc, L, Uij, Wij);
                U += Uij;
                W += Wij;
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
    assert(not x.empty());
    const uint32_t N(x.size() / 3);
    one_energy(x, ofs, rc, Urc, L, dU, dW);
    dU = -dU - (2.0 * static_cast<double>(N) - 1.0) * ULRC0;
    dW = -dW - (2.0 * static_cast<double>(N) - 1.0) * WLRC0;
}

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

uint32_t shuffle(vector<double>& x, double& U, double& W,
        const double kT, const double dxmax, const double L,
        const double rc, const double Urc)
{
    // shuffle atoms, MC
    uint32_t _3N = x.size();
    uint32_t N = _3N / 3;
    uint32_t Naccept(0);

    double Uold, Wold, Unew, Wnew;
    double dU, dW, beta_dU;
    vector<double> oldx(3);
    for (uint32_t ofs(0); ofs < _3N; ofs += 3) {
        one_energy(x, ofs, rc, Urc, L, Uold, Wold);
        for (int k(0); k < 3; ++k) {
            oldx[k] = x[ofs + k];
            x[ofs + k] += randomer::rand(-dxmax, dxmax);
            x[ofs + k] -= L * round(x[ofs + k] / L); // periodic condition
        }
        one_energy(x, ofs, rc, Urc, L, Unew, Wnew);

        dU = Unew - Uold;
        dW = Wnew - Wold;
        beta_dU = dU / kT;

        if (decide(beta_dU)) {
            U += dU;
            W += dW;
            Naccept += 1;
        }
        else {
            copy(oldx.begin(), oldx.end(), x.begin() + ofs);
        }
    } 
    return Naccept;
}

bool create(vector<double>& x, double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    uint32_t N(_3N / 3);
    double dU, dW;
    double dCB;
    vector<double> newx(randomer::vrand(3, -0.5 * L, 0.5 * L));
    potin(newx, x, rc, Urc, L, ULRC0, WLRC0, dU, dW);

    dCB = dU / kT - mu / kT  - log(V / (N + 1));
    if (decide(dCB)) {
        x.resize(_3N + 3);
        copy(newx.begin(), newx.end(), x.end() - 3);
        U += dU;
        W += dW;
        return true;
    }
    return false;
}

bool destruct(vector<double>& x, double& U, double& W,
        const double rc, const double Urc,
        const double L, const double V,
        const double kT, const double mu,
        const double ULRC0, const double WLRC0)
{
    uint32_t _3N(x.size());
    uint32_t N(_3N / 3);
    uint32_t ofs(randomer::choice(N) * 3);
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

void muVTMC() 
{
    // output config
    ioer::output_t out("");
    out.set_precision(4);

    // long range correction
    double sr3, sr9, ULRC0(0.0), WLRC0(0.0);
    sr3 = pow(1.0 / para.rc, 3);
    sr9 = sr3 * sr3 * sr3;
    /*
    ULRC0 = 8.0 / 9.0 * M_PI / para.V * (sr9 - 3.0 * sr3);
    WLRC0 = 32.0 / 9.0 * M_PI / para.V * (sr9 - 1.5 * sr3);
    */

    // shifted energy
    double Urc(0.0), Wrc(0.0);
    raw_pair_energy(para.rc * para.rc, Urc, Wrc);
    //Urc = 0.0;


    // prepare/load init configuration
    vector<double> x;
    if (para.prepinit) {
        x = randomer::vrand(3 * para.N0, -0.5 * para.L, 0.5 * para.L);
    }
    else {
        read_conf(x, para.conffile);
    }
    out.info("# init N = ", x.size() / 3);

    // calc init U & W
    double U, W;
    all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, U, W);
    out.info("# init U = ", U, " init W = ", W);

    // main MC part
    uint32_t N;
    uint32_t Naccept(0), Nmove(0);
    double Usum(0.0), Psum(0.0), rhosum(0.0);
    uint32_t Nsamp(0);

    // loop over steps
    double Unow, Pnow, rhonow;
    out.info("# start looping ... ");
    out.tabout("# Nsamp", "rho", "U/N", "P", "<rho>", "<U/N>", "<P>");
    for (uint32_t istep(0); istep < para.Nstep; ++istep) {

        // shuffle
        Naccept += shuffle(x, U, W, para.kT, para.dxmax, para.L, para.rc, Urc);
        N = x.size() / 3;
        Nmove += N;
        // create / destruct
        if (randomer::rand() < 0.5) {
            create(x, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
        }
        else {
            destruct(x, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
        }

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
    out.tabout(Nsamp, rhonow, Unow, Pnow);
    write_conf(x, para.conffile);

    //final output
    out.drawline('\n', 1);
    out.tabout("# <rho> = ", rhosum / Nsamp);
    out.tabout("# <U/N> = ", Usum / Nsamp);
    out.tabout("# <p> = ", Psum / Nsamp);
    out.tabout("# acc ratio = ", static_cast<double>(Naccept) / Nmove);
    out.drawline('\n', 1);

} // main


int main(int argc, char** argv) {

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
