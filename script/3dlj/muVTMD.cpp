#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
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
    /*
    x = randomer::vrand(3 * para.N0, -0.5 * para.L, 0.5 * para.L);
    */

    // init v
    v = randomer::maxwell_dist(mass, kT, N0);
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
        ("# Nstep", para.Nstep)
        ("# Anastep", para.Anastep)
        ("# mass", para.mass)
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
    double ULRC0, WLRC0, Urc;
    tail_correction(para.rc, para.V, para.LJmodel, Urc, ULRC0, WLRC0);
    out.info("# tail correction: Urc = ", Urc, " ULRC0 = ", ULRC0, " WLRC0 = ", WLRC0);

    // init configuration
    vector<double> x, v;
    double U, W;
    if (para.prepinit) {
        init_conf(x, v, para.mass, para.kT, para.L, para.N0);
    }
    else {
        read_conf(x, v, para.conffile);
    }
    all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, U, W, nullptr);
    out.info("# init configuration: N = ", x.size() / 3, " init U = ", U, " init W = ", W);

    // main MD part
    vector<double> vrec;
    double randnum;
    const double shuffle_frac(para.move_frac), create_frac(1.0 - (1.0 - shuffle_frac) * 0.5);
    uint64_t N;
    double Usum(0.0), Wsum(0.0), rhosum(0.0), Ksum(0.0), kT_sum(0.0), Psum(0.0);
    uint64_t Nsamp(0);
    out.info("# start MD looping ... ");
    out.info("# shuffle frac = ", shuffle_frac, " create frac = ", create_frac);
    out.tabout("# Nsamp", "<rho>", "<U>", "<P>", "<kT>", "rho", "U", "W");
    for (uint64_t istep(0); istep < para.Nstep; ++istep) {
        N = x.size() / 3;
        // evolve
        if (not x.empty()) {
            evolve(x, v, U, W, para.L, para.dt, para.mass, para.kT, para.nu,
                    para.rc, Urc, ULRC0, WLRC0);
            andersen_thermostat(v, para.mass, para.kT, para.nu, para.dt);
        }

        // exchange
        if (istep % para.K == 0) {
            if (randomer::rand() < 0.5) {
                // create
                vector<double> newx(randomer::vrand(3, -0.5 * para.Lc, 0.5 * para.Lc));
                vector<double> newv(randomer::maxwell_dist(para.mass, para.kT));
                create(x, v, newx, newv, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
            }
            else {
                // destruct
                vector<uint64_t> Vc_idx(get_Vc_idx(x, para.Lc));
                if (not Vc_idx.empty()) {
                    uint64_t ofs(randomer::choice(Vc_idx) * 3);
                    destruct(x, v, ofs, U, W, para.rc, Urc, para.L, para.V, para.kT, para.mu, ULRC0, WLRC0);
                }
            }
        }

        // statistics
        Nsamp += 1;
        N = x.size() / 3;
        if (N > 0) {
            //vrec.insert(vrec.end(), v.begin(), v.end());
            rhosum += N / para.V;
            Usum += U;
            Wsum += W;
            kT_sum += 2 * cal_Ek(v, para.mass) / N / 3.0;
        }
        // output & save
        if (istep % para.Anastep == 0) {
            // print statistics & save configure
            out.tabout(Nsamp, 
                    rhosum / Nsamp,
                    Usum / Nsamp / para.V / (rhosum / Nsamp),
                    rhosum / Nsamp * para.kT + Wsum / Nsamp / para.V,
                    kT_sum / Nsamp,
                    N / para.V,
                    U / N,
                    W / N
                    );
            write_conf(x, v, para.conffile);
        }
    } 
    out.tabout(Nsamp, 
            rhosum / Nsamp,
            Usum / Nsamp / para.V / (rhosum / Nsamp),
            rhosum / Nsamp * para.kT + Wsum / Nsamp / para.V,
            kT_sum / Nsamp,
            N / para.V,
            U / N,
            W / N
            );
    write_conf(x, v, para.conffile);
    /*
    // save vrec
    ioer::h5file_t f("vrec.dat", ios::out);
    f.create_dataset("vrec", vrec);
    f.close();
    */

    //final output
    const double avgrho(rhosum / Nsamp);
    const double avgU(Usum / Nsamp);
    const double avgW(Wsum / Nsamp);

    out.drawline('\n', 1);
    out.tabout("# kT = ", para.kT);
    out.tabout("# mu = ", para.mu);
    out.tabout("# <rho> = ", avgrho);
    out.tabout("# <U> = ", avgU / para.V / avgrho);
    out.tabout("# <P> = ", avgrho * para.kT + avgW / para.V);
    out.drawline('\n', 1);

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
