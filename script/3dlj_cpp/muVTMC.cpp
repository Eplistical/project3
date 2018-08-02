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
using namespace std;

// configuration I/O functions

void read_conf(vector<double>& x, const string& conffile) 
{
    ioer::h5file_t in(conffile, ios::in);
    uint32_t N;
    string program;
    in.read_attr("para", "N", N, "program", program);
    in.read_dataset("x", x);
    assert(N * 3 == x.size());
    assert(program == "muVTMC");
    in.close();
}

inline void write_conf(const vector<double>& x, const string& conffile) 
{
    ioer::h5file_t out(conffile, ios::out);
    uint32_t N(x.size() / 3);
    out.create_dataset("para", vector<double>{0.0});
    out.create_attr("para", "N", N, "program", "muVTMC");
    out.create_dataset("x", x);
    out.close();
}

// muVTMC

void muVTMC() 
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
        ("# rc", para.rc)
        ("# kT", para.kT)
        ("# mu", para.mu)
        ("# Nstep", para.Nstep)
        ("# dxmax", para.dxmax)
        ("# prepinit", para.prepinit)
        ("# N0", para.N0)
        ("# move_frac", para.move_frac)
        ("# conffile", para.conffile)
        ("# Anastep", para.Anastep)
        ("# random_seed", para.random_seed)
        ;
    out.info("# --- PROGRAM BEGINS --- ");
    out.info("# exe name: muVTMC ");

    // tail corrections
    double ULRC0, WLRC0, Urc;
    tail_correction(para.rc, para.V, para.LJmodel, Urc, ULRC0, WLRC0);
    out.info("# tail correction: Urc = ", Urc, " ULRC0 = ", ULRC0, " WLRC0 = ", WLRC0);

    // init configuration
    vector<double> x;
    double U, W;
    if (para.prepinit) {
        x = randomer::vrand(3 * para.N0, -0.5 * para.L, 0.5 * para.L);
    }
    else {
        read_conf(x, para.conffile);
    }
    all_energy(x, para.rc, Urc, para.L, ULRC0, WLRC0, U, W);
    out.info("# init configuration: N = ", x.size() / 3, " init U = ", U, " init W = ", W);

    // main MC part
    double randnum;
    const double shuffle_frac(para.move_frac), create_frac(1.0 - (1.0 - shuffle_frac) * 0.5);
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
        muVTMC();
        ioer::info("# ", cpu_timer.format(4));
    }
    return 0;
}
