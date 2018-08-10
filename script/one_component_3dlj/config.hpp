#ifndef _CONFIG_HPP
#define _CONFIG_HPP
#include <string>
#include <cassert>

#include "boost/program_options.hpp"
namespace po = boost::program_options;

// config
struct Para {
    // basic
    std::string LJmodel = "cs";
    double V = 512;
    double L = pow(V, 1.0 / 3.0); 
    double Vc = 1.0;
    double Lc;
    double rc = 2.5;

    double kT = 1.0; 
    double mu = -3.2;
    double mass = 1.0;

    // propagation
    uint64_t Nstep = 2e6;
    uint64_t Anastep = 1e4;

    // configuration
    bool prepinit = false;
    uint64_t N0 = 108;
    std::string conffile = "conf.dat";

    // MC related
    double dxmax = 1.0;
    double move_frac = 0.75;

    // MD related
    double dt = 0.005;
    double nu = 10;
    uint64_t K = 5;

    // others 
    uint64_t random_seed = 0;

    public:
    void show(ioer::output_t& out) {
        out.info("# --- CONFIG PARAMETERS --- ");
        out.keyval()
            ("# LJmodel", LJmodel)
            ("# V", V)
            ("# L", L)
            ("# Vc", Vc)
            ("# Lc", Lc)
            ("# rc", rc)
            ("# kT", kT)
            ("# mu", mu)
            ("# mass", mass)
            ("# Nstep", Nstep)
            ("# Anastep", Anastep)
            ("# prepinit", prepinit)
            ("# N0", N0)
            ("# conffile", conffile)
            ("# dxmax", dxmax)
            ("# move_frac", move_frac)
            ("# dt", dt)
            ("# nu", nu)
            ("# K", K)
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
        ("LJmodel", po::value<std::string>(&para.LJmodel), "specific LJ model version, c (cutoff) or cs (cutoff + shifted)")
        ("V", po::value<double>(&para.V), "box volume")
        ("Vc", po::value<double>(&para.Vc), "control volume, as a multiple of V")
        ("rc", po::value<double>(&para.rc), "cutoff distance for LJ")
        ("kT", po::value<double>(&para.kT), "temperature")
        ("mu", po::value<double>(&para.mu), "chemical potential")
        ("mass", po::value<double>(&para.mass), "mass")
        ("Nstep", po::value<uint64_t>(&para.Nstep), "total time step, non-negative integer")
        ("Anastep", po::value<uint64_t>(&para.Anastep), "analysis interval, non-negative integear")
        ("prepinit", po::value<bool>(&para.prepinit), "bool, if true, prepare initial configuration")
        ("N0", po::value<uint64_t>(&para.N0), "initial # of atoms to prepare, non-negative integer")
        ("conffile", po::value<std::string>(&para.conffile), "file for configuration")
        ("dxmax", po::value<double>(&para.dxmax), "MC max displacement on each direction")
        ("move_frac", po::value<double>(&para.move_frac), "MC probability to call move rather than create/destruct")
        ("dt", po::value<double>(&para.dt), "MD time step")
        ("nu", po::value<double>(&para.nu), "MD Anderen thermostat collision frequency")
        ("K", po::value<uint64_t>(&para.K), "MD interval to call create/remove, non-negative integer")
        ("random_seed", po::value<uint64_t>(&para.random_seed), "random seed, non-negative integer")
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
    assert(para.Vc <= 1.0 and para.Vc > 0.0);
    para.Vc *= para.V;
    para.Lc = pow(para.Vc, 1.0 / 3.0);

    if (vm.count("LJmodel")) {
        assert(para.LJmodel == "c" or para.LJmodel == "cs");
    }

    return true;
}


#endif 
