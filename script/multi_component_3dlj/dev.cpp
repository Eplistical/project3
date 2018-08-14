#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <map>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/crasher.hpp"
#include "misc/fmtstring.hpp"
#include "misc/ioer.hpp"
#include "boost/timer/timer.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "mcmd.hpp"
/*
#include "thermostat.hpp"
*/
using namespace std;

const string program_name("dev");

// configuration functions
void init_conf(vector<double>& x, vector<double>& v, 
        const uint64_t Ntype, vector<uint64_t>& type,
        const vector<double>& mass, const double kT,
        const vector<double>& L, const vector<uint64_t>& N0)
{
    // initialize N0 particles in the box with length L
    const uint64_t Ntot(sum(N0));

    // init type
    type.clear();
    for (uint64_t i(0); i < Ntype; ++i) {
        vector<uint64_t> tmp(N0[i], i);
        type.insert(type.end(), tmp.begin(), tmp.end());
    }
    //shuffle(type.begin(), type.end(), randomer::rng);

    // init x
    x.resize(Ntot * 3);
    uint64_t N_per_row(1);
    while (N_per_row * N_per_row * N_per_row < Ntot) {
        N_per_row += 1;
    }
    vector<double> spacing(L / N_per_row);
    uint64_t idx(0);
    for (uint64_t ix(0); ix < N_per_row; ++ix) {
        for (uint64_t iy(0); iy < N_per_row; ++iy) {
            for (uint64_t iz(0); iz < N_per_row; ++iz) {
                if (idx >= Ntot) {
                    continue;
                }
                else {
                    x[0 + idx * 3] = (ix + 0.5) * spacing[0] - 0.5 * L[0];
                    x[1 + idx * 3] = (iy + 0.5) * spacing[1] - 0.5 * L[0];
                    x[2 + idx * 3] = (iz + 0.5) * spacing[2] - 0.5 * L[0];
                    idx += 1;
                }
            }
        }
    }

    // init v
    v.clear();
    /*
    v.resize(Ntot * 3);
    for (uint64_t i(0); i < Ntot; ++i) {
        vector<double> tmp(randomer::maxwell_dist(mass[type[i]], kT));
        copy(tmp.begin(), tmp.begin() + 3, v.begin() + i * 3);
    }
    */
}

void read_conf(vector<double>& x, vector<uint64_t>& type, const string& conffile) 
{
    ioer::h5file_t in(conffile, ios::in);
    string program;
    in.read_attr("para", "program", program);
    in.read_dataset("x", x, "type", type);
    misc::crasher::confirm<>(program == program_name, "incompatibale program name!");
    in.close();
}

inline void write_conf(const vector<double>& x, const vector<uint64_t>& type, const string& conffile) 
{
    ioer::h5file_t out(conffile, ios::out);
    out.create_dataset("para", vector<double>{0.0});
    out.create_attr("para", "program", program_name);
    out.create_dataset("x", x, "type", type);
    out.close();
}

// run

void run() 
{
    // output config
    ioer::output_t out("");
    out.set_precision(6);

    // show para
    out.info("# --- CONFIG PARAMETERS --- ");
    out.keyval()
        ("# Ntype", para.Ntype)
        ("# LJmodel", para.LJmodel)
        ("# L", para.L)
        ("# V", para.V)
        ("# Lc", para.Lc)
        ("# Vc", para.Vc)
        ("# rc", para.rc)
        ("# sigma", para.sigma)
        ("# epsilon", para.epsilon)
        ("# kT", para.kT)
        ("# mu", para.mu)
        ("# mass", para.mass)
        ("# Nstep", para.Nstep)
        ("# Anastep", para.Anastep)
        ("# dxmax", para.dxmax)
        ("# prepinit", para.prepinit)
        ("# move_frac", para.move_frac)
        ("# conffile", para.conffile)
        ("# random_seed", para.random_seed)
        ;
    out.info("# --- PROGRAM BEGINS --- ");
    out.info("# program name: ", program_name);

    // tail corrections
    vector<double> Urc(para.Ntype * para.Ntype);
    vector<double> ULRC(para.Ntype * para.Ntype);
    vector<double> WLRC(para.Ntype * para.Ntype);
    tail_correction(para.Ntype, para.rc, para.sigma, para.epsilon, 
            para.V, para.LJmodel, Urc, ULRC, WLRC);
    out.info("# tail correction:");
    out.info("# Urc = ", Urc);
    out.info("# ULRC = ", ULRC);
    out.info("# WLRC = ", WLRC);

    // setup
    vector<double> x, v;
    vector<uint64_t> type;
    double U, W;
    if (para.prepinit) {
        init_conf(x, v, para.Ntype, type, para.mass, para.kT, para.L, para.N0);
    }
    else {
        read_conf(x, type, para.conffile);
    }
    all_energy(x, para.Ntype, type, para.sigma, para.epsilon, 
            para.rc, Urc, para.L, ULRC, WLRC, U, W, nullptr);

    vector<uint64_t> N(get_N(para.Ntype, type));
    uint64_t Ntot(sum(N));

    out.info("# init configuration: N = ", N);
    out.info("# init U = ", U, " init W = ", W);

    // observable map
    uint64_t Nsamp(0);
    map<string, double> obs;
    vector<string> obs_keys { "Usum", "Wsum"};
    for (int i(0); i < para.Ntype; ++i) {
        obs_keys.push_back(misc::fmtstring("rhosum%d", i));
    }
    for (const auto& key : obs_keys) {
        ioer::info(key);
        obs[key] = 0.0;
    }

    uint64_t Naccept(0), Nmove(0);
    out.info("# start MC looping ... ");
    out.tabout_nonewline("# Nsamp", "<U>", "<P>");
    for (int i(0); i < para.Ntype; ++i) {
        out.tabout_nonewline(misc::fmtstring("<rho%d>", i));
    }
    out.newline();

    double randnum;
    uint64_t Nc;
    N = get_N(para.Ntype, type);
    Ntot = sum(N);
    for (uint64_t istep(0); istep < para.Nstep; ++istep) {
        randnum = randomer::rand();

        // shuffle
        /*
        if (not x.empty()) {
            uint64_t idx(randomer::choice(Ntot));
            Naccept += shuffle(x, idx, para.Ntype, type, para.kT, para.dxmax, para.sigma, para.epsilon, 
                    para.rc, Urc, para.L, ULRC, WLRC, U, W);
        }
        */
        if (randnum < 0.5) {
            // create
            uint64_t newtype;
            vector<double> newx, newv;
            rand_in_Vc(para.Ntype, para.mass, para.kT, para.Lc, newtype, newx, newv);
            Nc = get_Nc(newtype, type, x, para.Lc);
            Naccept += create(x, v, type, newx, newv, newtype, para.Ntype, N, 
                    para.sigma, para.epsilon, para.rc, Urc, 
                    para.L, para.kT, para.mu, para.Vc, Nc,
                    ULRC, WLRC, U, W);
            Nc = get_Nc(newtype, type, x, para.Lc);
        }
        else {
            // destruct
            vector<uint64_t> Vc_idx(get_Vc_idx(x, para.Lc));
            if (not Vc_idx.empty()) {
                uint64_t idx(randomer::choice(Vc_idx));
                Nc = get_Nc(type[idx], type, x, para.Lc);
                Naccept += destruct(x, v, type, idx, 
                        para.Ntype, N,
                        para.sigma, para.epsilon, para.rc, Urc,
                        para.L, para.kT, para.mu, para.Vc, Nc,
                        ULRC, WLRC, U, W);
                Nc = get_Nc(type[idx], type, x, para.Lc);
            }
        }
        Nmove += 1;

        // statistics
        Nsamp += 1;
        N = get_N(para.Ntype, type);
        Ntot = sum(N);
        if (Ntot > 0) {
            obs["Usum"] += U;
            obs["Wsum"] += W;
            for (int i(0); i < para.Ntype; ++i) {
                obs[misc::fmtstring("rhosum%d", i)] += N[i] / para.V;
            }
        }
        // output & save
        if (istep % para.Anastep == 0) {
            // get rho sum
            double rhosum(0.0);
            for (int i(0); i < para.Ntype; ++i) {
                rhosum += obs[misc::fmtstring("rhosum%d", i)];
            }

            // output
            out.tabout_nonewline(Nsamp, 
                    obs["Usum"] / para.V / rhosum,
                    rhosum / Nsamp * para.kT + obs["Wsum"] / Nsamp / para.V
                    );

            for (int i(0); i < para.Ntype; ++i) {
                out.tabout_nonewline(obs[misc::fmtstring("rhosum%d", i)] / Nsamp);
            }
            out.newline();

            // save conf
            write_conf(x, type, para.conffile);
        }
    }
    // get rho sum
    double rhosum(0.0);
    for (int i(0); i < para.Ntype; ++i) {
        rhosum += obs[misc::fmtstring("rhosum%d", i)];
    }

    // output
    out.tabout_nonewline(Nsamp, 
            obs["Usum"] / para.V / rhosum,
            rhosum / Nsamp * para.kT + obs["Wsum"] / Nsamp / para.V
            );

    for (int i(0); i < para.Ntype; ++i) {
        out.tabout_nonewline(obs[misc::fmtstring("rhosum%d", i)] / Nsamp);
    }
    out.newline();

    // save conf
    write_conf(x, type, para.conffile);

    out.tabout("# acc ratio = ", static_cast<double>(Naccept) / Nmove);
} 

// main

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    else {
        randomer::seed(para.random_seed);

        boost::timer::cpu_timer cpu_timer;
        run();
        ioer::info("# ", cpu_timer.format(4));
    }
    return 0;
}
