#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <map>
#include "misc/vector.hpp"
#include "misc/randomer.hpp"
#include "misc/ioer.hpp"
#include "boost/timer/timer.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "mcmd.hpp"
#include "thermostat.hpp"
using namespace std;

const string program_name("muVTMD");

// configuration functions
void init_conf(vector<double>& x, vector<double>& v,
        const double mass, const double kT,
        const double L, const uint64_t N0) 
{
    // initialize N0 particles in the box with length L
    const uint64_t _3N(N0 * 3);

    // init x
    x.resize(_3N);
    uint64_t N_per_row(1);
    while (N_per_row * N_per_row * N_per_row < N0) {
        N_per_row += 1;
    }
    double spacing(L / N_per_row);
    uint64_t idx(0);
    for (int ix(0); ix < N_per_row; ++ix) {
        for (int iy(0); iy < N_per_row; ++iy) {
            for (int iz(0); iz < N_per_row; ++iz) {
                if (idx >= N0) {
                    continue;
                }
                else {
                    x[0 + idx * 3] = (ix + 0.5) * spacing - 0.5 * L;
                    x[1 + idx * 3] = (iy + 0.5) * spacing - 0.5 * L;
                    x[2 + idx * 3] = (iz + 0.5) * spacing - 0.5 * L;
                    idx += 1;
                }
            }
        }
    }
    // init v
    //v = randomer::maxwell_dist(mass, kT, N0);
    v.resize(N0 * 3);
    for (uint64_t i(0); i < N0; ++i) {
        vector<double> tmp(randomer::maxwell_dist(mass, kT));
        copy(tmp.begin(), tmp.begin() + 3, v.begin() + i * 3);
    }
}

void read_conf(vector<double>& x, vector<double>& v, const string& conffile) 
{
    ioer::h5file_t in(conffile, ios::in);
    uint64_t N;
    string program;
    in.read_attr("para", "N", N, "program", program);
    if (N > 0) {
        in.read_dataset("x", x);
        in.read_dataset("v", v);
    }
    else {
        x.resize(0);
        v.resize(0);
    }
    assert(N * 3 == x.size());
    assert(N * 3 == v.size());
    assert(program == program_name);
    in.close();
}

inline void write_conf(const vector<double>& x, const vector<double>& v, const string& conffile) 
{
    ioer::h5file_t out(conffile, ios::out);
    uint64_t N(x.size() / 3);
    out.create_dataset("para", vector<double>{0.0});
    out.create_attr("para", "N", N, "program", program_name);
    out.create_dataset("x", x);
    out.create_dataset("v", v);
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
        ("# LJmodel", para.LJmodel)
        ("# V", para.V)
        ("# L", para.L)
        ("# Vc", para.Vc)
        ("# Lc", para.Lc)
        ("# rc", para.rc)
        ("# kT", para.kT)
        ("# mu", para.mu)
        ("# mass", para.mass)
        ("# Nstep", para.Nstep)
        ("# Anastep", para.Anastep)
        ("# dt", para.dt)
        ("# nu", para.nu)
        ("# K", para.K)
        ("# prepinit", para.prepinit)
        ("# conffile", para.conffile)
        ("# random_seed", para.random_seed)
        ;
    out.info("# --- PROGRAM BEGINS --- ");
    out.info("# program name: ", program_name);

    // tail corrections
    double ULRC, WLRC, Urc;
    tail_correction(para.rc, para.V, para.LJmodel, Urc, ULRC, WLRC);
    out.info("# tail correction:");
    out.info("# Urc = ", Urc);
    out.info("# ULRC = ", ULRC);
    out.info("# WLRC = ", WLRC);

    // init configuration
    vector<double> x, v;
    double U, W;
    if (para.prepinit) {
        init_conf(x, v, para.mass, para.kT, para.L, para.N0);
    }
    else {
        read_conf(x, v, para.conffile);
    }
    all_energy(x, para.rc, Urc, para.L, ULRC, WLRC, U, W, nullptr);
    out.info("# init configuration: N = ", x.size() / 3);
    out.info("# init U = ", U);
    out.info("# init W = ", W);

    // main MD part
    double randnum;
    uint64_t N;

    // observable map
    uint64_t Nsamp(0);
    map<string, double> obs;
    vector<string> obs_keys { "rhosum", "Usum", "Wsum", "kTsum"};
    for (const auto& key : obs_keys) {
        obs[key] = 0.0;
    }

    // loop
    out.info("# start MD looping ... ");
    out.tabout("# Nsamp", "<U>", "<P>", "<kT>", "<rho>");
    for (uint64_t istep(0); istep < para.Nstep; ++istep) {
        N = x.size() / 3;
        // evolve
        if (not x.empty()) {
            evolve(x, v, U, W, para.L, para.dt, para.mass, para.kT, para.nu,
                    para.rc, Urc, ULRC, WLRC);
            andersen_thermostat(v, para.mass, para.kT, para.nu, para.dt);
        }

        // exchange
        if (istep % para.K == 0) {
            vector<uint64_t> Vc_idx(get_Vc_idx(x, para.Lc));
            uint64_t Nc(Vc_idx.size());

            if (randomer::rand() < 0.5) {
                // create
                vector<double> newx(randomer::vrand(3, -0.5 * para.Lc, 0.5 * para.Lc));
                vector<double> newv(randomer::maxwell_dist(para.mass, para.kT));
                create(x, v, newx, newv, U, W, para.rc, Urc, para.L, para.Vc, Nc, para.kT, para.mu, ULRC, WLRC);
            }
            else {
                // destruct
                if (not Vc_idx.empty()) {
                    uint64_t ofs(randomer::choice(Vc_idx) * 3);
                    destruct(x, v, ofs, U, W, para.rc, Urc, para.L, para.Vc, Nc, para.kT, para.mu, ULRC, WLRC);
                }
            }
        }

        // statistics
        Nsamp += 1;
        N = x.size() / 3;
        if (N > 0) {
            obs["rhosum"] += N / para.V;
            obs["Usum"] += U;
            obs["Wsum"] += W;
            obs["kTsum"] += 2.0 * cal_Ek(v, para.mass) / 3.0 / N;
        }
        // output & save
        if (istep % para.Anastep == 0) {

            // output
            out.tabout(Nsamp,
                    obs["Usum"] / para.V / obs["rhosum"],
                    (obs["rhosum"] * obs["kTsum"] / Nsamp + obs["Wsum"] / para.V) / Nsamp,
                    obs["kTsum"] / Nsamp,
                    obs["rhosum"] / Nsamp
                    );
            write_conf(x, v, para.conffile);
        }
    } 
    out.tabout(Nsamp,
            obs["Usum"] / para.V / obs["rhosum"],
            obs["rhosum"] / Nsamp * obs["kTsum"] / Nsamp + obs["Wsum"] / Nsamp / para.V,
            obs["kTsum"] / Nsamp,
            obs["rhosum"] / Nsamp
            );
    write_conf(x, v, para.conffile);

    //final output
    out.info("# --- PROGRAM ENDS --- ");
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
